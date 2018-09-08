using ArgParse
using CSV
using DataFrames
using Mamba
using MicrobiomeMixedModels
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
            help = "YAML defining the nodes to monitor through MCMC."
        "--inits"
            help = "YAML defining initial values of stochastic nodes."
        "--hyper"
            help = "YAML defining hyperparameters of stochastic node priors."
        "--factors"
            arg_type = Int
            help = "Number of factors to use in fitting the MIMIX model."
        "--no-factors"
            help = "Use the MIMIX w/o Factors model (ignores --factors)."
            action = :store_true
        "--permanova"
            help = "Run PERMANOVA with vegan::adonis() in R."
            action = :store_true
        "--post-pred-check"
            help = "Perform posterior predictive checks, requires at least one of :Y, :ϕ, or :θ are monitored."
    end
    return parse_args(s)
end

function read_data(dir; L=0)
    println("Reading X.csv")
    X = readcsv(joinpath(dir, "X.csv"), Int)
    println("Reading Y.csv")
    Y = readcsv(joinpath(dir, "Y.csv"), Int)
    Y = Y[:, 1:100]
    println(size(Y))
    println("Reading Z.csv")
    Z = readcsv(joinpath(dir, "Z.csv"), Int)
    N, K = size(Y)
    m = sum(Y, 2)
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
        bf = find(any(Z .== level, 1))
        @assert length(bf) == 1
        data[:blocking_factor][level] = bf[1]
    end

    return data
end

args = parse_commandline()

input = abspath(args["input"])
@assert isdir(input)

output = abspath(args["output"])

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
    if args["no-factors"] | factors == 0
        model_type = MIMIXNoFactors()
        factors = 0
    elseif factors > 0
        model_type = MIMIX(factors)
    else
        ValueError("--factors requires positive integer or --no-factors flag must be given")
    end

    monitor_conf = load_config(abspath(args["monitor"]))
    hyper_conf = load_config(abspath(args["hyper"]))
    inits_conf = load_config(abspath(args["inits"]))

    model = get_model(model_type, monitor_conf, hyper_conf)

    data = read_data(input; L = factors)

    inits = get_inits(model_type, inits_conf, data)
    inits = [inits for _ in 1:args["chains"]]

    println("Beginning MCMC sampling")
    mcmc_kwargs = Dict(Symbol(key) => args[key] for key in ["burnin", "thin", "chains"])
    sim = mcmc(model, data, inits, args["iters"]; mcmc_kwargs...)

    posterior_out = joinpath(output, "posterior-samples.jls")
    println("Saving posterior samples to $posterior_out")
    write(posterior_out, sim)
end