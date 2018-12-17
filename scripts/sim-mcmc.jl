push!(LOAD_PATH, "/Users/neal/Projects/mimix/src")

using ArgParse
using CSV
using DataFrames
using Mamba
using MicrobiomeMixedModels
using Random
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
        "--seed"
            arg_type = Int
            default = 1
            help = "Reseed the random number generator."
        "--data"
            nargs = '*'
            help = "YAML defining the settings for artificial data generation."
        "--monitor"
            nargs = '*'
            help = "YAML defining the nodes to monitor through MCMC."
        "--inits"
            nargs = '*'
            help = "YAML defining initial values of stochastic nodes."
        "--hyper"
            nargs = '*'
            help = "YAML defining hyperparameters of stochastic node priors."
        "--factors"
            arg_type = Int
            default = 0
            help = "Number of factors to use in fitting the MIMIX model."
        "--loadings"
            default = "DL"
            help = "Distributions on loadings (G for Gaussian priors, DL for Dirichlet-Laplace priors)."
        "--permanova"
            help = "Run PERMANOVA with vegan::adonis() in R."
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
    μ::Vector=collect(range(1, stop=-1, length=K)), block_cor::Float64=0.9, error_cor::Float64=0.9,
    form::String="grouped", dense::Float64=0.1, block_var::Float64=1.0, error_var::Float64=1.0,
    a_support::Vector{Float64}=collect(0.01:0.01:0.5), seed::Int=1)

    @assert length(μ) == K
    @assert iseven(N)
    @assert N % q == 0

    Random.seed!(seed)

    # Fixed effects due to single treatment
    X = transpose([ones(div(N, 2))... zeros(div(N, 2))...])
    β = zeros(K)
    if form == "grouped"
        ngroup = round(Int64, dense / 0.05)
        group = shuffle(collect(1:5:K))[1:ngroup]
        for g in group
            v = rand(Uniform(-2, 2))
            v += 1.0sign(v)
            for i in g:min(g + 4, K)
                β[i] = v
            end
        end
    elseif form == "random"
        nindividual = round(Int64, dense * K)
        individual = shuffle(collect(1:K))[1:nindividual]
        for i in individual
            v = rand(Uniform(-2, 2))
            v += 1.0sign(v)
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
    ϕ = e ./ sum(e, dims=2)
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
        bf = findall(vec(any(Z .== level, dims=1)))
        @assert length(bf) == 1
        data[:blocking_factor][level] = bf[1]
    end

    return data, truth
end

args = parse_commandline()

output = abspath(args["output"])
if !isdir(output)
    mkdir(output)
end

@assert 0 <= args["seed"]   "Seed must be non-negative"

if args["permanova"]
    data_conf = load_config(AbstractString[abspath(data_path) for data_path in args["data"]])

    seed = args["seed"]
    data, truth = generate_data(; seed = seed, data_conf...)

    println("Beginning PERMANOVA test with vegan::adonis in R.")
    sim = permanova(data)
    global_results = DataFrame(
        reject_global_null = sim[:pval] < 0.05,
        dense = data_conf[:dense],
        block_var = data_conf[:block_var],
        error_var = data_conf[:error_var],
        form = data_conf[:form]
    )

    global_test_path = joinpath(output, "global-test.tsv")
    println("Writing global test results to $global_test_path")
    CSV.write(global_test_path, global_results, delim='\t')
else  # mimix
    @assert 0 < args["iters"]   "Iters must be positive"
    @assert 0 <= args["burnin"] "Burn-in must be non-negative"
    @assert 0 < args["thin"]    "Thin must be positive"
    @assert 0 < args["chains"]  "Chains must be positive"

    factors = args["factors"]
    loadings = args["loadings"]
    if (factors > 0) & (loadings == "DL")
        model_type = MIMIX(factors)
    elseif (factors > 0) & (loadings == "G")
        model_type = MIMIXGaussian(factors)
    elseif factors == 0
        model_type = MIMIXNoFactors()
    else
        ValueError("--factors requires a non-negative integer.")
    end

    monitor_conf = load_config(AbstractString[abspath(data_path) for data_path in args["monitor"]])
    hyper_conf = load_config(AbstractString[abspath(data_path) for data_path in args["hyper"]])
    data_conf = load_config(AbstractString[abspath(data_path) for data_path in args["data"]])
    inits_conf = load_config(AbstractString[abspath(data_path) for data_path in args["inits"]])

    model = get_model(model_type, monitor_conf, hyper_conf)

    seed = args["seed"]
    data, truth = generate_data(; L = factors, seed = seed, data_conf...)

    inits = get_inits(model_type, inits_conf, data)
    inits = [inits for _ in 1:args["chains"]]

    println("Beginning MCMC simulation")
    Random.seed!(args["seed"])
    mcmc_kwargs = Dict(Symbol(key) => args[key] for key in ["burnin", "thin", "chains"])
    sim = mcmc(model, data, inits, args["iters"]; mcmc_kwargs...)

    # summarize global test results
    sim_omega = sim[:, [startswith(name, 'ω') for name in sim.names], :].value
    num_included_each_iter = sum(sim_omega, dims=2)
    post_prob_inclusion = mean(num_included_each_iter .> 0.0)
    global_results = DataFrame(
        reject_global_null = post_prob_inclusion > 0.9,
        dense = data_conf[:dense],
        block_var = data_conf[:block_var],
        error_var = data_conf[:error_var],
        form = data_conf[:form]
    )

    global_test_path = joinpath(output, "global-test.tsv")
    println("Writing global test results to $global_test_path")
    CSV.write(global_test_path, global_results, delim='\t')

    # summarize local parameter estimates
    sim_beta_names = sim.names[[startswith(name, 'β') for name in sim.names]]
    sim_beta = sim[:, sim_beta_names, :]
    results = DataFrame(mamba_name = sim_beta_names)
    nodes = Symbol[]
    vals = Float64[]
    for name in results[:mamba_name]
        for (node, value) in truth
            if startswith(name, String(node))
                push!(nodes, node)
                if '[' in name
                    index = name[collect(findfirst("[", name))[1]:end]
                    index = strip(index, ['[', ']'])
                    index = parse.(Int, split(index, ','))
                    push!(vals, value[index...])
                else
                    push!(vals, value)
                end
            end
        end
    end
    results[:mamba_node] = nodes
    results[:value] = vals

    post_summary = summarystats(sim_beta)
    post_quantiles = quantile(sim_beta)
    results[:mean] = post_summary.value[:, 1]
    for (i, q) in enumerate(post_quantiles.colnames)
        results[Symbol(q)] = post_quantiles.value[:, i]
    end

    local_estimates_path = joinpath(output, "local-estimates.tsv")
    println("Writing local parameter estimates to $local_estimates_path")
    CSV.write(local_estimates_path, results, delim='\t')
end