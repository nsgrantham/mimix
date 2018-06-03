using ArgParse
using CSV
using DataFrames
using Mamba
using MicrobiomeMixedModels
using YAML

function parse_commandline()
    s = ArgParseSettings("Fit a model to artificially-generated data.")
    @add_arg_table s begin
        "output"
            help = "File to which simulation results are written."
        "--iters", "-i"
            arg_type = Int
            default = 500
            help = "Number of MCMC iterations."
        "--burnin", "-b"
            arg_type = Int
            default = 0
            help = "Number of initial MCMC iterations to remove as burn-in."
        "--thin", "-t"
            arg_type = Int
            default = 1
            help = "Retain one of every 'thin' iterations."
        "--chains", "-c"
            arg_type = Int
            default = 1
            help = "Number of MCMC chains to run."
        "--data"
            help = "YAML defining the settings for artificial data generation."
        "--monitor"
            help = "YAML defining the nodes to monitor through MCMC."
        "--inits"
            help = "YAML defining initial values of stochastic nodes."
        "--hyper"
            help = "YAML defining hyperparameters of stochastic node priors."
        "--factors"
            arg_type = Int
            help = "Number of factors to use in fitting the MIMIX model."
        "--no-factors"
            help = "Ignore --factors and use the MIMIX w/o Factors model"
            action = :store_true
    end
    return parse_args(s)
end

function ar1_covar(d::Int, cor::Float64, var::Float64)
    @assert var > zero(var)
    @assert abs(cor) <= one(cor)
    Σ = zeros(d, d)
    for j in 1:d
        for i in 1:d
            Σ[i, j] = var * cor^abs(i - j)
        end
    end
    Σ
end

function generate_data(;
    N::Int=40, K::Int=100, L::Int=-1, q::Int=5, m_min::Int=2500, m_max::Int=5000,
    μ::Vector=collect(linspace(1, -1, K)), block_cor::Float64=0.9, error_cor::Float64=0.9,
    dense::Float64=0.1, block_var::Float64=1.0, error_var::Float64=1.0,
    a_support::Vector{Float64}=collect(0.01:0.01:0.5))

    @assert length(μ) == K
    @assert iseven(N)
    @assert N % q == 0

    # Fixed effects due to single treatment
    X = transpose([ones(div(N, 2))... zeros(div(N, 2))...])
    β = zeros(K)
    ngroup = round(Int64, dense / 0.05)
    group = shuffle(collect(1:5:K))[1:ngroup]
    for g in group
        v = rand(Uniform(-2, 2))
        v += 1.0sign(v)
        for i in g:min(g + 4, K)
            β[i] = v
        end
    end

    # Random effects due to block
    Z = reshape(repeat(collect(1:q), outer=[div(N, q)]), N, 1)
    γ = rand(MvNormal(ar1_covar(K, block_cor, block_var)), q)

    # Unconstrained transformed probabilities
    Σ = ar1_covar(K, error_cor, error_var)
    θ = zeros(N, K)
    for i in 1:N
        θ[i, :] = rand(MvNormal(μ + β * X[i] + γ[:, Z[i]], Σ))
    end

    # Multivariate count data
    e = exp.(θ)
    ϕ = e ./ sum(e, 2)
    m = sample(m_min:m_max, N)
    Y = zeros(Int, (N, K))
    for i in 1:N
        Y[i, :] = rand(Multinomial(m[i], ϕ[i, :]))
    end


    # Define and return dicts of data (for mcmc) and truth (for validation)
    truth = Dict{Symbol, Any}(
        :β => β
    )

    data = Dict{Symbol, Any}(
        :X => X,
        :Y => Y,
        :Z => Z,
        :N => N,
        :K => K,
        :m => m,
        :q => q,
        :L => L
    )
    data[:p] = size(X, 2)
    data[:num_blocking_factors] = size(Z, 2)
    data[:blocking_factor] = Dict{Int, Int}()
    for level in unique(Z)
        bf = find(any(Z .== level, 1))
        @assert length(bf) == 1
        data[:blocking_factor][level] = bf[1]
    end

    return data, truth
end

args = parse_commandline()

@assert 0 < args["iters"]   "Iters must be positive"
@assert 0 <= args["burnin"] "Burn-in must be non-negative"
@assert 0 < args["thin"]    "Thin must be positive"
@assert 0 < args["chains"]  "Chains must be positive"

factors = args["factors"]
if args["no-factors"]
    mimix = MIMIXNoFactors()
    factors = -1
elseif factors > 0
    mimix = MIMIX(factors)
else
    ValueError("--factors requires positive integer or --no-factors flag must be given")
end

monitor_conf = load_config(abspath(args["monitor"]))
hyper_conf = load_config(abspath(args["hyper"]))
data_conf = load_config(abspath(args["data"]))
inits_conf = load_config(abspath(args["inits"]))

model = get_model(mimix, monitor_conf, hyper_conf)
data, truth = generate_data(; L = factors, data_conf...)
#data, truth = generate_data(; L = factors)
inits = get_inits(mimix, inits_conf, data)
inits = [inits for _ in 1:args["chains"]]

mcmc_kwargs = Dict(Symbol(key) => args[key] for key in ["burnin", "thin", "chains"])
sim = mcmc(model, data, inits, args["iters"]; mcmc_kwargs...)

# summarize simulation results in DataFrame
results = DataFrame(MambaName = sim.names)
nodes = Symbol[]
values = Float64[]
for name in results[:MambaName]
    for (node, value) in truth
        if startswith(name, String(node))
            push!(nodes, node)
            if '[' in name
                index = name[search(name, '['):end]
                index = strip(index, ['[', ']'])
                index = parse.(split(index, ','))
                push!(values, value[index...])
            else
                push!(values, value)
            end
        end
    end
end
results[:MambaNode] = nodes
results[:Value] = values

post_summary = summarystats(sim)
post_quantiles = quantile(sim)
results[:Mean] = post_summary.value[:, 1]
for (i, q) in enumerate(post_quantiles.colnames)
    results[Symbol(q)] = post_quantiles.value[:, i]
end

print(results)
#writetable(abspath(args["output"]), results)
CSV.write(abspath(args["output"]), results, delim='\t')
