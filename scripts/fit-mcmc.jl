push!(LOAD_PATH, "/Users/neal/Projects/mimix/src")

using ArgParse
using CSV
using DataFrames
using DelimitedFiles
using LinearAlgebra
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

    if args["post-pred-check"] && !monitor_conf[:θ] && (!monitor_conf[:θ_mean] || !monitor_conf[:θ_var])
        error("Posterior prediction checks require monitoring θ or both θ_mean and θ_var, please edit the monitor config.")
    end

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
    monitor_θ_mean = pop!(monitor_conf, :θ_mean)
    monitor_θ_var = monitor_conf[:θ_var]
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
        # ΛtΛs = Matrix{Float64}[]
        # for i in 1:size(Λ_values, 1)
        #     Λ_mat = reshape(Λ_values[i, :], data[:L], data[:K])
        #     push!(ΛtΛs, transpose(Λ_mat) * Λ_mat)
        # end
        # ΛtΛ_postmean = sum(ΛtΛs, dims=3)
        # ΛtΛ_postmean ./= size(Λ_values, 1)
        # post_path = joinpath(output, "LambdatLambda-postmean.tsv")
        # println("Saving posterior mean of Lambda' * Lambda to $post_path")
        # writedlm(post_path, ΛtΛ_postmean)
    end

    if args["post-pred-check"]
        # Posterior predictive checks on sparsity, overdispersion, and ecological diversity (alpha and beta)

        function shannon_diversity(x::Vector{T}) where T <: Real
            p = x ./ sum(x)
            p = p[p .> 0.]
            return -sum(p .* log.(p))
        end
        
        function simpson_diversity(x::Vector{T}) where T <: Real
            p = x ./ sum(x)
            return sum(abs2, p)
        end

        function braycurtis_diversity(x::Vector{T}, y::Vector{T}) where T <: Real
            c = zero(eltype(x))
            s = zero(eltype(x))
            for i in eachindex(x, y)
                c += min(x[i], y[i])
                s += x[i] + y[i]
            end
            return 2c / s
        end

        function jaccard_diversity(x::Vector{T}, y::Vector{T}) where T <: Real
            p = 0
            q = 0
            r = 0
            z = zero(T)
            for i in eachindex(x, y)
                if (x[i] > z) && (y[i] > z)
                    p += 1
                elseif (x[i] > z) && (y[i] == z)
                    q += 1
                elseif (x[i] == z) && (y[i] > z)
                    r += 1
                end
            end
            return p / (p + q + r)
        end


        function calculate_mean_beta_div_by_group(Y, groups, diversity_fn)
            @assert size(Y, 1) == length(groups)
            mean_div = zeros(maximum(group_ids))
            for group_id in unique(group_ids)
                Y_group = Y[group_id .== group_ids, :]
                pairwise_div = Float64[]
                for i in 1:size(Y_group, 1)
                    for j in i+1:size(Y_group, 1)
                        push!(pairwise_div, diversity_fn(Y_group[i, :], Y_group[j, :]))
                    end
                end
                mean_div[group_id] = mean(pairwise_div)
            end
            return mean_div
        end

        matchrows(A::Matrix, B::Matrix) = mapslices(findall, hcat([all(A .== transpose(B[i, :]), dims=2) for i in 1:size(B, 1)]...), dims=2) 

        full_combos = hcat(data[:X], data[:Z])
        group_ids = vec(matchrows(full_combos, unique(full_combos, dims=1)))  # each group is unique combo of treatment, site, and blocks
        
        obs_max_count = vec(maximum(data[:Y], dims=2) ./ sum(data[:Y], dims=2))
        obs_mean_count = vec(var(data[:Y], dims=2))
        obs_prop_eq_zero = vec(mean(data[:Y] .== 0.0, dims=2))
        obs_prop_leq_one = vec(mean(data[:Y] .<= 1.0, dims=2))
        obs_prop_leq_two = vec(mean(data[:Y] .<= 2.0, dims=2))
        obs_shannon_div = mapslices(shannon_diversity, data[:Y], dims=2)
        obs_simpson_div = mapslices(simpson_diversity, data[:Y], dims=2)
        obs_braycurtis_div = reshape(calculate_mean_beta_div_by_group(data[:Y], group_ids, braycurtis_diversity), maximum(group_ids), 1)
        obs_jaccard_div = reshape(calculate_mean_beta_div_by_group(data[:Y], group_ids, jaccard_diversity), maximum(group_ids), 1)
        writedlm(joinpath(output, "obs-max-count.tsv"), obs_max_count)
        writedlm(joinpath(output, "obs-mean-count.tsv"), obs_mean_count)
        writedlm(joinpath(output, "obs-prop-eq-zero.tsv"), obs_prop_eq_zero)
        writedlm(joinpath(output, "obs-prop-leq-one.tsv"), obs_prop_leq_one)
        writedlm(joinpath(output, "obs-prop-leq-two.tsv"), obs_prop_leq_two)
        writedlm(joinpath(output, "obs-shannon-div.tsv"), obs_shannon_div)
        writedlm(joinpath(output, "obs-simpson-div.tsv"), obs_simpson_div)
        writedlm(joinpath(output, "obs-braycurtis-div.tsv"), obs_braycurtis_div)
        writedlm(joinpath(output, "obs-jaccard-div.tsv"), obs_jaccard_div)

        post_pred_checks = Dict{String, Matrix}()

        if monitor_θ
            θ_post = convert(Matrix, get_post(sim, data, :θ))
            post_pred_checks["theta_mcmc"] = θ_post
        end
        
        if monitor_θ_mean && monitor_θ_var
            θ_mean_post = convert(Matrix, get_post(sim, data, :θ_mean))
            θ_var_post = convert(Matrix, get_post(sim, data, :θ_var))
            iters = size(θ_mean_post, 1) 
            θ_post = zeros(iters, data[:N] * data[:K])
            for iter in 1:iters
                θ_mean = reshape(θ_mean_post[iter, :], data[:N], data[:K])
                θ_new = [rand(MvNormal(θ_mean[i, :], sqrt.(θ_var_post[iter, :]))) for i in 1:data[:N]]
                θ_post[iter, :] = vec(transpose(hcat(θ_new...)))
            end
            post_pred_checks["theta_prior"] = θ_post
        end

        for post_pred_check in post_pred_checks
            name, θ_post = post_pred_check

            n_iter = size(θ_post, 1)
            max_count = zeros(n_iter, data[:N])
            mean_count = zeros(n_iter, data[:N])
            prop_eq_zero = zeros(n_iter, data[:N])
            prop_leq_one = zeros(n_iter, data[:N])
            prop_leq_two = zeros(n_iter, data[:N])
            shannon_div = zeros(n_iter, data[:N])
            simpson_div = zeros(n_iter, data[:N])
            braycurtis_div = zeros(n_iter, maximum(group_ids))
            jaccard_div = zeros(n_iter, maximum(group_ids))
            for iter in 1:n_iter
                ϕ = clr(reshape(θ_post[iter, :], data[:N], data[:K]))
                Y_pred_iter = transpose(hcat([rand(Multinomial(data[:m][i], ϕ[i, :])) for i in 1:data[:N]]...))
                max_count[iter, :] = vec(maximum(Y_pred_iter, dims=2) ./ sum(Y_pred_iter, dims=2))
                mean_count[iter, :] = vec(var(Y_pred_iter, dims=2))
                prop_eq_zero[iter, :] = vec(mean(Y_pred_iter .== 0.0, dims=2))
                prop_leq_one[iter, :] = vec(mean(Y_pred_iter .<= 1.0, dims=2))
                prop_leq_two[iter, :] = vec(mean(Y_pred_iter .<= 2.0, dims=2))
                shannon_div[iter, :] = mapslices(shannon_diversity, Y_pred_iter, dims=2)
                simpson_div[iter, :] = mapslices(simpson_diversity, Y_pred_iter, dims=2)
                braycurtis_div[iter, :] = calculate_mean_beta_div_by_group(Y_pred_iter, group_ids, braycurtis_diversity) 
                jaccard_div[iter, :] = calculate_mean_beta_div_by_group(Y_pred_iter, group_ids, jaccard_diversity) 
            end

            if !isdir(joinpath(output, name))
                mkdir(joinpath(output, name))
            end

            writedlm(joinpath(output, name, "post-pred-max-count.tsv"), max_count)
            writedlm(joinpath(output, name, "post-pred-mean-count.tsv"), mean_count)
            writedlm(joinpath(output, name, "post-pred-prop-eq-zero.tsv"), prop_eq_zero)
            writedlm(joinpath(output, name, "post-pred-prop-leq-one.tsv"), prop_leq_one)
            writedlm(joinpath(output, name, "post-pred-prop-leq-two.tsv"), prop_leq_two)
            writedlm(joinpath(output, name, "post-pred-shannon-div.tsv"), shannon_div)
            writedlm(joinpath(output, name, "post-pred-simpson-div.tsv"), simpson_div)
            writedlm(joinpath(output, name, "post-pred-braycurtis-div.tsv"), braycurtis_div)
            writedlm(joinpath(output, name, "post-pred-jaccard-div.tsv"), jaccard_div)
        end
    end
end