###  Cross-validation comparing MIMIX & MIMIX w/ Factors

# Set up the validation in three parts:
#   1. Data (the data to allocate to different folds)
#   2. Partition (a function to partition the data into different folds)
#   3. Metrics (summary metrics to compare performance of different models),
#   4. Methods (models to evaluate over different folds and calculate metrics)


### Data

@everywhere data = Dict{Symbol, Any}(
    :Y => readcsv("data/Y.csv", Int),
    :X => readcsv("data/X.csv", Int),
    :Z => readcsv("data/Z.csv", Int)
)


### Partition

@everywhere using StatsBase

@everywhere function partition(data, folds::Int)
    Y = data[:Y]
    n, K = size(Y)
    m = vec(sum(Y, 2))
    Y_train = [zeros(n, K) for j in 1:folds]
    Y_test = [zeros(n, K) for j in 1:folds]
    for i in 1:n
        y = Y[i, :]
        taxa = zeros(Int, m[i])
        fold = zeros(Int, m[i])
        obs_idx = 1
        for (k, yk) in enumerate(y)
            for obs in 1:yk
                taxa[obs_idx] = k
                fold[obs_idx] = sample(1:folds)
                obs_idx += 1
            end
        end
        for j in 1:folds
            test_set = fold .== j
            Y_test[j][i, :] = fit(Histogram, taxa[test_set], 0:K).weights
            Y_train[j][i, :] = fit(Histogram, taxa[!test_set], 0:K).weights
        end
    end
    [(Y_train[j], Y_test[j], data[:X], data[:Z]) for j in 1:folds]
end


### Metrics

@everywhere using DataStructures: OrderedDict

metrics = OrderedDict{Symbol, String}()
for i in 1:size(data[:Y], 1)
    key = Symbol(@sprintf("ll%03d", i))
    metrics[key] = "Log-likelihood evaluated at observation $i"
end


### Methods

@everywhere include("models.jl")
@everywhere settings = Dict{Symbol, Any}(
    :monitor => [:ϕ],
    :iters => 20000,
    :burnin => 10000,
    :thin => 2,
    :epsilon => 0.05,
    :steps => 16
)

@everywhere function summarize(mc::Mamba.ModelChains, Y_test)
    # Estimate occurrence probability vector ϕ
    n, K = size(Y_test)
    ϕ_post = mc[:, vec(["ϕ[$i,$k]" for i in 1:n, k in 1:K]), :].value
    ϕ_post = vcat([ϕ_post[:, :, i] for i in 1:size(ϕ_post, 3)]...)
    ϕ_postmean = reshape(vec(mean(ϕ_post, 1)), n, K)
    ϕ_postmean ./= sum(ϕ_postmean, 2)

    # Evaluate log-likelihood for each observation i = 1,...,n
    m_test = vec(sum(Y_test, 2))
    results = Dict{Symbol, Any}()
    for i in 1:n
        key = Symbol(@sprintf("ll%03d", i))
        results[key] = logpdf(Multinomial(m_test[i], ϕ_postmean[i, :]), Y_test[i, :])
    end
    results
end

@everywhere methods = Dict{Symbol, Function}(
    :MIMIX => (Y_train, Y_test, X, Z) ->
        summarize(fit(MIMIX(166), Y_train, X, Z; settings...), Y_test),
    :MIMIXNoFactors => (Y_train, Y_test, X, Z) ->
        summarize(fit(MIMIXNoFactors(), Y_train, X, Z; settings...), Y_test)
)


### Run the cross-validation!

@everywhere include("utils.jl")
validate(partition, data, metrics, methods;
         dir="results", filename="validate.csv", folds=5)
