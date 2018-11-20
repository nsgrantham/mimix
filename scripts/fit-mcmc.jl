push!(LOAD_PATH, "/Users/neal/Projects/mimix/src")

using ArgParse
using CSV
using DataFrames
using DelimitedFiles
using Mamba
using MicrobiomeMixedModels
using Random
using YAML

function parse_commandline()
    s = ArgParseSettings("Fit a model to real-world data.")
    @add_arg_table s begin
        "input"
            help = "Directory from which to read X.csv, Y.csv, and Z.csv."
        "output"
            help = "Directory to which results are written."
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
        "--post-pred-check"
            help = "Perform posterior predictive checks, requires at least one of :Y, :ϕ, or :θ are monitored."
            action = :store_true
    end
    return parse_args(s)
end

function read_data(dir; L=0)
    println("Reading X.csv")
    X = convert(Matrix{Int}, CSV.read(joinpath(dir, "X.csv"), datarow=1))
    println("Reading Y.csv")
    Y = convert(Matrix{Int}, CSV.read(joinpath(dir, "Y.csv"), datarow=1))
    println("Reading Z.csv")
    Z = convert(Matrix{Int}, CSV.read(joinpath(dir, "Z.csv"), datarow=1))
    N, K = size(Y)
    m = sum(Y, dims=2)
    q = maximum(vec(Z))
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

    return data
end

args = parse_commandline()

input = abspath(args["input"])
@assert isdir(input)

output = abspath(args["output"])

Random.seed!(args["seed"])

if args["permanova"]
    data = read_data(input)

    println("Beginning PERMANOVA test with vegan::adonis in R.")
    sim = permanova(data)
    global_results = DataFrame(pval = sim[:pval])

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
        ValueError("--factors requires a non-negative integer")
    end

    monitor_conf = load_config(AbstractString[abspath(data_path) for data_path in args["monitor"]])
    hyper_conf = load_config(AbstractString[abspath(data_path) for data_path in args["hyper"]])
    inits_conf = load_config(AbstractString[abspath(data_path) for data_path in args["inits"]])

    model = get_model(model_type, monitor_conf, hyper_conf)

    data = read_data(input; L = factors)

    inits = get_inits(model_type, inits_conf, data)
    inits = [inits for _ in 1:args["chains"]]

    println("Beginning MCMC sampling")
    mcmc_kwargs = Dict(Symbol(key) => args[key] for key in ["burnin", "thin", "chains"])
    sim = mcmc(model, data, inits, args["iters"]; mcmc_kwargs...)

    monitor_Y = pop!(monitor_conf, :Y)
    monitor_ϕ = pop!(monitor_conf, :ϕ)
    monitor_θ = pop!(monitor_conf, :θ)
    monitor_Λ = pop!(monitor_conf, :Λ)

    for (param, monitor) in monitor_conf
        if monitor
            post = get_post(sim, data, param)
            param = String(param)
            for (letter, name) in MicrobiomeMixedModels.latin
                if startswith(param, letter)
                    param = replace(param, letter => name)
                end
            end
            post_path = joinpath(output, "$param.tsv")
            println("Saving posterior samples of $param to $post_path")
            CSV.write(post_path, post, delim='\t')
        end
    end

    if monitor_Λ
        # Lambda is large, so rather than save every sample we save the posterior mean
        Λ = get_post(sim, data, :Λ)
        Λ_values = convert(Matrix, Λ)
        Λ_postmean = reshape(vec(mean(Λ_values, dims=1)), data[:L], data[:K])
        post_path = joinpath(output, "Lambda-postmean.tsv")
        println("Saving posterior mean of Lambda to $post_path")
        writedlm(post_path, Λ_postmean)

        # We can also store the posterior mean of Lambda' * Lambda
        ΛtΛs = Matrix{Float64}[]
        for i in 1:size(Λ_values, 1)
            Λ_mat = reshape(Λ_values[i, :], data[:L], data[:K])
            push!(ΛtΛs, transpose(Λ_mat) * Λ_mat)
        end
        ΛtΛ_postmean = sum(ΛtΛs, dims=3)
        ΛtΛ_postmean ./= size(Λ_values, 1)
        post_path = joinpath(output, "LambdatLambda-postmean.tsv")
        println("Saving posterior mean of Lambda' * Lambda to $post_path")
        writedlm(post_path, ΛtΛ_postmean)
    end

    if args["post-pred-check"]
        # Posterior predictive checks
        obs_max_count = vec(maximum(data[:Y], dims=2) ./ sum(data[:Y], dims=2))
        obs_prop_eq_zero = vec(mean(data[:Y] .== 0.0, dims=2))
        obs_prop_leq_one = vec(mean(data[:Y] .<= 1.0, dims=2))
        obs_prop_leq_two = vec(mean(data[:Y] .<= 2.0, dims=2))
        writedlm(joinpath(output, "obs-max-count.tsv"), obs_max_count)
        writedlm(joinpath(output, "obs-prop-eq-zero.tsv"), obs_prop_eq_zero)
        writedlm(joinpath(output, "obs-prop-leq-one.tsv"), obs_prop_leq_one)
        writedlm(joinpath(output, "obs-prop-leq-two.tsv"), obs_prop_leq_two)

        if monitor_θ
            θ_post = convert(Matrix, get_post(sim, data, :θ))
        end

        ϕ = zeros(data[:N], data[:K])
        Y_pred_iter = zeros(data[:N], data[:K])
        n_iter = size(θ_post, 1)
        max_count = zeros(n_iter, data[:N])
        prop_eq_zero = zeros(n_iter, data[:N])
        prop_leq_one = zeros(n_iter, data[:N])
        prop_leq_two = zeros(n_iter, data[:N])
        for i in 1:n_iter
            ϕ[:, :] = clr(reshape(θ_post[i, :], data[:N], data[:K]))
            Y_pred_iter[:, :] = hcat([rand(Multinomial(data[:m][i], ϕ[i, :])) for i in 1:data[:N]]...)'
            max_count[i, :] = vec(maximum(Y_pred_iter, dims=2) ./ sum(Y_pred_iter, dims=2))
            prop_eq_zero[i, :] = vec(mean(Y_pred_iter .== 0.0, dims=2))
            prop_leq_one[i, :] = vec(mean(Y_pred_iter .<= 1.0, dims=2))
            prop_leq_two[i, :] = vec(mean(Y_pred_iter .<= 2.0, dims=2))
        end

        writedlm(joinpath(output, "post-pred-max-count.tsv"), max_count)
        writedlm(joinpath(output, "post-pred-prop-eq-zero.tsv"), prop_eq_zero)
        writedlm(joinpath(output, "post-pred-prop-leq-one.tsv"), prop_leq_one)
        writedlm(joinpath(output, "post-pred-prop-leq-two.tsv"), prop_leq_two)
    end
end