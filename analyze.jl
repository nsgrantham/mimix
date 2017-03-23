include("models.jl")

Y = readcsv("data/Y.csv", Int)
X = readcsv("data/X.csv", Int)
Z = readcsv("data/Z.csv", Int)
n, K = size(Y)
p = size(X, 2)
num_blocking_factors = size(Z, 2)

## MCMC settings for MIMIX & MIMIX w/o Factors
settings = Dict{Any, Any}(
    :iters => 20000,
    :burnin => 10000,
    :thin => 2,
    :hmc_verbose => false,  # if true, will print HMC accept rate to console
    :epsilon => 0.05,       # step size, tune for accept b/w 25-45%
    :steps => 16            # typically powers of 2, tune for accept b/w 25-45%
)


### MIMIX

srand(881)
L = n

sim = fit(MIMIX(L), Y, X, Z; settings...,
          monitor=[:β, :b_var, :g_var, :Λ, :θ_var, :ω, :θ])

## Save posterior samples for params of interest
mimix_dir = joinpath("results", "analyze", "mimix")
function writepost(filename, names)
    post = sim[:, names, :].value
    post = vcat([post[:, :, i] for i in 1:size(post, 3)]...)
    writecsv(joinpath(mimix_dir, filename), post)
end
writepost("b-var.csv", ["b_var[$j]" for j in 1:p])
writepost("g-var.csv", ["g_var[$r]" for r in 1:num_blocking_factors])
writepost("theta-var.csv", ["θ_var[$k]" for k in 1:K])
writepost("omega1.csv", ["ω[$l,1]" for l in 1:L])
writepost("omega2.csv", ["ω[$l,2]" for l in 1:L])
writepost("omega3.csv", ["ω[$l,3]" for l in 1:L])
writepost("beta1.csv", ["β[$k,1]" for k in 1:K])
writepost("beta2.csv", ["β[$k,2]" for k in 1:K])
writepost("beta3.csv", ["β[$k,3]" for k in 1:K])

# Lambda is large, so rather than save every sample we save the posterior mean
Λ = sim[:, vec(["Λ[$l,$k]" for l in 1:L, k in 1:K]), :].value
Λ = vcat([Λ[:, :, i] for i in 1:size(Λ, 3)]...)
Λ_postmean = reshape(vec(mean(Λ, 1)), L, K)
writecsv(joinpath(mimix_dir, "Lambda-postmean.csv"), Λ_postmean)

# We can also store Lambda' * Lambda
ΛtΛ_postmean = zeros(K, K)
Λ_mat = zeros(L, K)
for i in 1:size(Λ, 1)
    Λ_mat[:, :] = reshape(Λ[i, :], L, K)
    ΛtΛ_postmean += Λ_mat' * Λ_mat
end
ΛtΛ_postmean ./= size(Λ, 1)
writecsv(joinpath(mimix_dir, "LambdatLambda-postmean.csv"), ΛtΛ_postmean)

# Posterior predictive checks
m = sum(Y, 2)
obs_max_count = vec(maximum(Y, 2) ./ sum(Y, 2))
obs_prop_eq_zero = vec(mean(Y .== 0.0, 2))
obs_prop_leq_one = vec(mean(Y .<= 1.0, 2))
obs_prop_leq_two = vec(mean(Y .<= 2.0, 2))
writecsv(joinpath(mimix_dir, "obs-max-count.csv"), obs_max_count)
writecsv(joinpath(mimix_dir, "obs-prop-eq-zero.csv"), obs_prop_eq_zero)
writecsv(joinpath(mimix_dir, "obs-prop-leq-one.csv"), obs_prop_leq_one)
writecsv(joinpath(mimix_dir, "obs-prop-leq-two.csv"), obs_prop_leq_two)

θ_post = sim[:, vec(["θ[$i,$k]" for i in 1:n, k in 1:K]), :].value
θ_post = vcat([θ_post[:, :, i] for i in 1:size(θ_post, 3)]...)
ϕ = zeros(n, K)
Y_pred_iter = zeros(n, K)
n_iter = size(θ_post, 1)
max_count = zeros(n_iter, n)
prop_eq_zero = zeros(n_iter, n)
prop_leq_one = zeros(n_iter, n)
prop_leq_two = zeros(n_iter, n)
for i in 1:n_iter
    ϕ[:, :] = clr(reshape(θ_post[i, :], n, K))
    Y_pred_iter[:, :] = hcat([rand(Multinomial(m[i], ϕ[i, :])) for i in 1:n]...)'
    max_count[i, :] = vec(maximum(Y_pred_iter, 2) ./ sum(Y_pred_iter, 2))
    prop_eq_zero[i, :] = vec(mean(Y_pred_iter .== 0.0, 2))
    prop_leq_one[i, :] = vec(mean(Y_pred_iter .<= 1.0, 2))
    prop_leq_two[i, :] = vec(mean(Y_pred_iter .<= 2.0, 2))
end

writecsv(joinpath(mimix_dir, "post-pred-max-count.csv"), max_count)
writecsv(joinpath(mimix_dir, "post-pred-prop-eq-zero.csv"), prop_eq_zero)
writecsv(joinpath(mimix_dir, "post-pred-prop-leq-one.csv"), prop_leq_one)
writecsv(joinpath(mimix_dir, "post-pred-prop-leq-two.csv"), prop_leq_two)
sim = 0

### MIMIX w/o Factors

srand(881)
sim = fit(MIMIXNoFactors(), Y, X, Z; settings...,
          monitor=[:β, :β_var, :γ_var, :θ_var, :ω, :θ])

## Save posterior samples for params of interest
mimix_dir = joinpath("results", "analyze", "mimix-no-factors")
function writepost(filename, names)
    post = sim[:, names, :].value
    post = vcat([post[:, :, i] for i in 1:size(post, 3)]...)
    writecsv(joinpath(mimix_dir, filename), post)
end
writepost("beta-var.csv", ["β_var[$j]" for j in 1:p])
writepost("gamma-var.csv", ["γ_var[$r]" for r in 1:num_blocking_factors])
writepost("theta-var.csv", ["θ_var[$k]" for k in 1:K])
writepost("omega1.csv", ["ω[$k,1]" for k in 1:K])
writepost("omega2.csv", ["ω[$k,2]" for k in 1:K])
writepost("omega3.csv", ["ω[$k,3]" for k in 1:K])
writepost("beta1.csv", ["β[$k,1]" for k in 1:K])
writepost("beta2.csv", ["β[$k,2]" for k in 1:K])
writepost("beta3.csv", ["β[$k,3]" for k in 1:K])

# Posterior predictive checks
m = sum(Y, 2)
obs_max_count = vec(maximum(Y, 2) ./ sum(Y, 2))
obs_prop_eq_zero = vec(mean(Y .== 0.0, 2))
obs_prop_leq_one = vec(mean(Y .<= 1.0, 2))
obs_prop_leq_two = vec(mean(Y .<= 2.0, 2))
writecsv(joinpath(mimix_dir, "obs-max-count.csv"), obs_max_count)
writecsv(joinpath(mimix_dir, "obs-prop-eq-zero.csv"), obs_prop_eq_zero)
writecsv(joinpath(mimix_dir, "obs-prop-leq-one.csv"), obs_prop_leq_one)
writecsv(joinpath(mimix_dir, "obs-prop-leq-two.csv"), obs_prop_leq_two)

θ_post = sim[:, vec(["θ[$i,$k]" for i in 1:n, k in 1:K]), :].value
θ_post = vcat([θ_post[:, :, i] for i in 1:size(θ_post, 3)]...)
ϕ = zeros(n, K)
Y_pred_iter = zeros(n, K)
n_iter = size(θ_post, 1)
max_count = zeros(n_iter, n)
prop_eq_zero = zeros(n_iter, n)
prop_leq_one = zeros(n_iter, n)
prop_leq_two = zeros(n_iter, n)
for i in 1:n_iter
    ϕ[:, :] = clr(reshape(θ_post[i, :], n, K))
    Y_pred_iter[:, :] = hcat([rand(Multinomial(m[i], ϕ[i, :])) for i in 1:n]...)'
    max_count[i, :] = vec(maximum(Y_pred_iter, 2) ./ sum(Y_pred_iter, 2))
    prop_eq_zero[i, :] = vec(mean(Y_pred_iter .== 0.0, 2))
    prop_leq_one[i, :] = vec(mean(Y_pred_iter .<= 1.0, 2))
    prop_leq_two[i, :] = vec(mean(Y_pred_iter .<= 2.0, 2))
end

writecsv(joinpath(mimix_dir, "post-pred-max-count.csv"), max_count)
writecsv(joinpath(mimix_dir, "post-pred-prop-eq-zero.csv"), prop_eq_zero)
writecsv(joinpath(mimix_dir, "post-pred-prop-leq-one.csv"), prop_leq_one)
writecsv(joinpath(mimix_dir, "post-pred-prop-leq-two.csv"), prop_leq_two)


### PERMANOVA w/ Bray-Curtis dissimilarity

using RCall
R"""
set.seed(881)
library(vegan)
z <- as.data.frame($Z)
z[] <- lapply(z, as.factor)
colnames(z) <- c("site", "block")
x <- as.data.frame($X)
x[] <- lapply(x, as.factor)
colnames(x) <- c("exclusion", "supplement", "combo")
y <- data.frame($Y)
dat <- cbind.data.frame(x, z)
sim <- adonis(y ~ site/block + supplement*exclusion, method="bray", data=dat, nperm=9999)
print(sim)
"""

