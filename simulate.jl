###  Simulation study comparing MIMIX, MIMIX w/o Factors, & PERMANOVA

# Set up the study in four parts:
#   1. Generate (a function to make data and true parameter values),
#   2. Factors (different features of the data gen process to vary),
#   3. Metrics (summary metrics to compare performance of different models),
#   4. Methods (models to evaluate at factor combinations and calculate metrics)


### Generate

@everywhere using Distributions

@everywhere function ar1(d::Int, cor::Float64, var::Float64)
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

@everywhere function generate(;
    N::Int=40, K::Int=100, q::Int=5, m_min::Int=2500, m_max::Int=5000, 
    μ::Vector=collect(linspace(1, -1, K)), block_cor::Float64=0.9, error_cor::Float64=0.9,
    dense::Float64=0.1, block_var::Float64=1.0, error_var::Float64=1.0)

    @assert length(μ) == K
    @assert iseven(N)
    @assert N % q == 0

    ## Fixed effects due to single treatment
    X = [ones(div(N, 2))... zeros(div(N, 2))...]'
    β = zeros(K)
    ngroup = round(Int64, dense / 0.05)
    group = shuffle(collect(1:5:K))[1:ngroup]
    for g in group
        v = rand(Uniform(-2, 2))
        v += 1.0sign(v)
        for i in g:min(g+4, K)
            β[i] = v
        end
    end

    ## Random effects due to block
    Z = reshape(repeat(collect(1:q), outer=[div(N, q)]), N, 1)
    γ = rand(MvNormal(ar1(K, block_cor, block_var)), q)

    ## Unconstrained transformed probabilities
    Σ = ar1(K, error_cor, error_var)
    θ = zeros(N, K)
    for i in 1:N
        θ[i, :] = rand(MvNormal(μ + β * X[i] + γ[:, Z[i]], Σ))
    end

    ## Multivariate count data
    e = exp(θ)
    ϕ = e ./ sum(e, 2)
    m = sample(m_min:m_max, N)
    Y = zeros(Int, (N, K))
    for i in 1:N
        Y[i, :] = rand(Multinomial(m[i], ϕ[i, :]))
    end
    
    (Y, X, Z, β)
end


###  Factors

factors = Dict{Symbol, Any}(
    :dense => [0.0, 0.05, 0.1],
    :error_var => [1.0, 4.0, 9.0],
    :block_var => [1.0, 4.0]
)

###  Metrics

@everywhere using DataStructures: OrderedDict

metrics = OrderedDict(
    :pval => "Global p-value (only from PERMANOVA)",
    :reject => "Reject global null?",
    :rmse => "Root mean squared error",
    :mae => "Mean absolute error",
    :coverage => "95% credible interval coverage",
    :tpr => "TPR",
    :fdr => "FDR",
    :tnr => "TNR",
    :f1 => "F1 score"
)


###  Methods

@everywhere include("models.jl")
@everywhere R"library(vegan)"
@everywhere settings = Dict{Symbol, Any}(
    :monitor => [:β, :ω],
    :iters => 10000,
    :burnin => 5000,
    :thin => 1,
    :epsilon => 0.15,
    :steps => 16
)

@everywhere function summarize(mc::Mamba.ModelChains, β)
    # Global test
    sim_omega = mc[:, [startswith(s, 'ω') for s in mc.names], :].value
    num_included_each_iter = sum(sim_omega, 2)
    post_prob_inclusion = mean(num_included_each_iter .> 0.0)
    reject = post_prob_inclusion > 0.9

    # Local tests & parameter estimation
    sim_beta = mc[:, ["β[$k,1]" for k in 1:length(β)], :]
    β_hat = Mamba.summarystats(sim_beta).value[:, 1]
    mae = mean(abs(β_hat - β))
    rmse = sqrt(mean((β_hat - β).^2))
    quant = Mamba.quantile(sim_beta).value
    β_lb = quant[:, 1]
    β_ub = quant[:, 5]
    coverage = mean(β_lb .<= β .<= β_ub)
    is_signif = β .!= 0.0
    is_predicted = !(β_lb .<= 0.0 .<= β_ub)
    tp = sum(is_predicted[is_signif])
    tn = sum(!is_predicted[!is_signif])
    fp = sum(is_predicted[!is_signif])
    fn = sum(!is_predicted[is_signif])
    tpr = tp / (tp + fn)
    fdr = fp / (tp + fp)
    tnr = tn / (tn + fp)
    f1  = 2tp / (2tp + fp + fn)

    Dict{Symbol, Any}(
        :rmse => rmse,
        :mae => mae,
        :coverage => coverage,
        :reject => reject,
        :tpr => tpr,
        :fdr => fdr,
        :tnr => tnr,
        :f1  => f1
    )
end


@everywhere methods = Dict{Symbol, Function}(
    :MIMIX => (Y, X, Z, β) -> 
        summarize(fit(MIMIX(40), Y, X, Z; settings...), β),
    :MIMIXNoFactors => (Y, X, Z, β) -> 
        summarize(fit(MIMIXNoFactors(), Y, X, Z; settings...), β),
    :PERMANOVA => (Y, X, Z, β) -> 
        begin
            R"""
            z <- data.frame($Z)
            z[] <- lapply(z, as.factor)
            x <- data.frame($X)
            x[] <- lapply(x, as.factor)
            sim <- adonis($Y ~ $Z + $X, method="bray", nperm=999)
            pval <- sim$aov.tab$`Pr(>F)`[1]
            """
            @rget pval
            Dict{Symbol, Any}(:pval => pval)
        end
)


### Run the simulation!

@everywhere include("utils.jl")
simulate(generate, factors, metrics, methods; 
         dir=joinpath("results", "simulate"), reps=50)
