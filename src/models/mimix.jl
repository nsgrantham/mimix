immutable MIMIX
    factors::Int
    function MIMIX(factors::Int)
        @assert factors > zero(factors)
        new(factors)
    end
end

function fit(mm::MIMIX, Y, X, Z;
             iters::Int=1000, monitor::Vector{Symbol}=Symbol[],
             epsilon::Float64=0.1, steps::Int=16, hmc_verbose::Bool=false,
             kwargs...)
    d = datadict(mm, Y, X, Z)
    i = inits(mm, d)
    m = model(mm, monitor)
    s = samplers(mm, epsilon, steps, hmc_verbose)
    setsamplers!(m, s)
    sim = mcmc(m, d, i, iters; kwargs...)
    return sim
end

function datadict(mm::MIMIX, Y, X, Z)
    d = Dict{Symbol, Any}(
        :Y => Y,
        :X => X,
        :Z => Z
    )
    d[:N], d[:K] = size(d[:Y])
    d[:m] = vec(sum(d[:Y], 2))
    d[:p] = size(d[:X], 2)
    d[:num_blocking_factor] = size(d[:Z], 2)
    d[:blocking_factor] = Dict{Int, Int}()
    for level in unique(d[:Z])
        bf = find(any(d[:Z] .== level, 1))
        @assert length(bf) == 1
        d[:blocking_factor][level] = bf[1]
    end
    d[:q] = maximum(d[:Z])
    d[:arr_a] = collect(0.01:0.01:0.5)
    d[:L] = mm.factors
    return d
end

function inits(::MIMIX, d::Dict{Symbol, Any})
    θ_init = iclr(proportionalize(d[:Y]))
    i = [
        Dict{Symbol, Any}(
            :Y => d[:Y],
            :θ => θ_init,
            :θ_var => vec(var(θ_init, 1)),
            :μ => vec(mean(θ_init, 1)),
            :μ_var => 1.0,
            :Λ => eye(d[:L], d[:K]),
            :F => zeros(d[:N], d[:L]),
            :ψ => ones(d[:L], d[:K]),
            :ξ => ones(d[:L], d[:K]),
            :τ => ones(d[:L]),
            :ν => 0.5,
            :b_full => zeros(d[:L], d[:p]),
            :g => zeros(d[:L], d[:q]),
            :ω => ones(d[:L], d[:p]),
            :π => [0.5 for j in 1:d[:p]],
            :b_var => ones(d[:p]),
            :g_var => ones(d[:num_blocking_factor]),
            :idx_a => ones(d[:L])
        )
    ]
    return i
end

function model(mm::MIMIX, monitor::Vector{Symbol}=Symbol[])
    Model(
        # High-dimensional counts

        Y = Stochastic(2,
            (m, ϕ, N) -> MultivariateDistribution[
                Multinomial(m[i], ϕ[i, :]) for i in 1:N
            ],
            :Y in monitor
        ),

        ϕ = Logical(2,
            (θ) -> clr(θ),
            :ϕ in monitor
        ),

        # Sample-specific random effects

        θ = Stochastic(2,
            (θ_mean, θ_var, N) -> MultivariateDistribution[
                MvNormal(θ_mean[i, :], sqrt(θ_var)) for i in 1:N
            ],
            :θ in monitor
        ),

        θ_mean = Logical(2,
            (μ, ΛF) -> ΛF .+ μ',
            :θ_mean in monitor
        ),

        θ_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            :θ_var in monitor
        ),

        # Population mean

        μ = Stochastic(1,
            (K, μ_var) -> MvNormal(K, sqrt(μ_var)),
            :μ in monitor
        ),

        μ_var = Stochastic(
            () -> InverseGamma(1, 1),
            :μ_var in monitor
        ),

        # Full-dimensional fixed effect estimates (monitor these for OTU-level inference!)

        β = Logical(2,
            (Λ, b) -> Λ' * b,
            :β in monitor
        ),

        # Factor analysis w/ loadings matrix Λ (L by K) and scores F (N by L)

        ΛF = Logical(2,
            (Λ, F) -> F * Λ,
            :ΛF in monitor
        ),

        Λ = Stochastic(2,
            (L, K, λ_var) -> MultivariateDistribution[
                MvNormal(zeros(K), sqrt(λ_var[l, :])) for l in 1:L
            ],
            :Λ in monitor
        ),

        λ_var = Logical(2,
            (ψ, τ, ξ) -> ψ .* (τ.^2 .* ξ.^2),
            :λ_var in monitor
        ),

        F = Stochastic(2,
            (F_mean, N) -> MultivariateDistribution[
                MvNormal(F_mean[i, :], 1.0) for i in 1:N
            ],
            :F in monitor
        ),

        F_mean = Logical(2,
            (Xb, sumZg) -> Xb + sumZg,
            :F_mean in monitor
        ),

        # Terms for Dirichlet-Laplace prior on factors

        ψ = Stochastic(2,
            (L, K) -> UnivariateDistribution[
                Exponential(2.0) for l in 1:L, k in 1:K
            ],
            :ψ in monitor
        ),

        ξ = Stochastic(2,
            (a, L, K) -> MultivariateDistribution[
                Dirichlet(fill(a[l], K)) for l in 1:L
            ],
            :ξ in monitor
        ),

        τ = Stochastic(1,
            (a, ν, L, K) -> UnivariateDistribution[
                Gamma(K * a[l], 1/ν) for l in 1:L
            ],
            :τ in monitor
        ),

        ν = Stochastic(
            () -> Gamma(1, 1),
            :ν in monitor
        ),

        a = Logical(1,
            (arr_a, idx_a, L) -> Float64[
                arr_a[convert(Int, idx_a[l])] for l in 1:L
            ],
            :a in monitor
        ),

        idx_a = Stochastic(1,
            (arr_a, L) -> UnivariateDistribution[
                Categorical(length(arr_a)) for l in 1:L
            ],
            :idx_a in monitor
        ),

        # Blocking factor random effects

        sumZg = Logical(2,
            (Zg) -> squeezesum(Zg, 3)',
            :sumZg in monitor
        ),

        Zg = Logical(3,
            (Z, g) -> g[:, Z],
            :Zg in monitor
        ),

        g = Stochastic(2,
            (q, L, g_var, blocking_factor) -> UnivariateDistribution[
                Normal(0.0, sqrt(g_var[blocking_factor[r]])) for l in 1:L, r in 1:q
            ],
            :g in monitor
        ),

        g_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            :g_var in monitor
        ),

        # Fixed treatment effects w/ spike-and-slab variable selection

        Xb = Logical(2,
            (X, b) -> X * b',
            :Xb in monitor
        ),

        b = Logical(2,
            (b_full, ω) -> ω .* b_full,
            :b in monitor
        ),

        b_full = Stochastic(2,
            (p, L, b_var) -> MultivariateDistribution[
                MvNormal(zeros(p), sqrt(b_var)) for l in 1:L
            ],
            :b_full in monitor
        ),

        b_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            :b_var in monitor
        ),

        ω = Stochastic(2,
            (p, L, π) -> UnivariateDistribution[
                Bernoulli(π[j]) for l in 1:L, j in 1:p
            ],
            :ω in monitor
        ),

        π = Stochastic(1,
            (L) -> Beta(1, L),
            :π in monitor
        )
    )
end

function samplers(::MIMIX, epsilon::Float64, steps::Int, hmc_verbose::Bool)

    HMC_θ = Sampler([:θ],
        (Y, m, θ, θ_mean, θ_var, N, K) ->
            begin
                acc = 0
                for i in 1:N
                    eps = rand(Exponential(epsilon))
                    yi = Y[i, :]
                    mi = m[i]
                    u1 = θ[i, :]
                    θi_mean = θ_mean[i, :]
                    logf0, grad0 = logf1, grad1 = logfgrad(u1, θi_mean, θ_var, yi, mi)
                    v0 = v1 = randn(length(u1))
                    v1 += 0.5eps * grad0
                    for step in 1:steps
                        u1 += eps * v1
                        logf1, grad1 = logfgrad(u1, θi_mean, θ_var, yi, mi)
                        v1 += eps * grad1
                    end
                    v1 -= 0.5eps * grad1
                    v1 *= -1.0
                    Kv0 = 0.5sumabs2(v0)
                    Kv1 = 0.5sumabs2(v1)
                    if rand() < exp((logf1 - Kv1) - (logf0 - Kv0))
                        acc += 1
                        θ[i, :] = u1
                    end
                end
                if hmc_verbose
                    println("HMC accept rate: $(round(100acc/N, 1))%, ")
                end
                θ
            end
    )

    Gibbs_μ = Sampler([:μ],
        (θ, θ_var, μ, μ_var, ΛF, N, K) ->
            begin
                Sigma = 1 ./ (N ./ θ_var .+ 1 / μ_var)
                mu = Sigma .* vec(sum(θ - ΛF, 1)) ./ θ_var
                for k in 1:K
                    μ[k] = rand(Normal(mu[k], sqrt(Sigma[k])))
                end
                μ
            end
    )

    Gibbs_Λ = Sampler([:Λ],
        (θ, θ_var, μ, Λ, λ_var, F, K, L) ->
            begin
                FtF = F' * F
                B = (F' * (θ .- μ')) ./ θ_var'
                Sigma = zeros(L, L)
                for k in 1:K
                    Sigma[:, :] = inv(cholfact(spdiagm(1 ./ λ_var[:, k]) + FtF ./ θ_var[k]))
                    Λ[:, k] = rand(MvNormal(Sigma * B[:, k], Hermitian(Sigma)))
                end
                Λ
            end
    )

    Gibbs_ψ = Sampler([:ψ],
        (Λ, ξ, ψ, τ, L, K) ->
            begin
                A = abs(Λ)
                boundbelow!(A)
                for k in 1:K
                    for l in 1:L
                        ψ[l, k] = 1 / rand(InverseGaussian(ξ[l, k] * τ[l] / A[l, k]))
                    end
                end
                ψ
            end
    )

    Gibbs_ξ = Sampler([:ξ],
        (Λ, a, ν, L, K) ->
            begin
                M = 2abs(Λ)
                boundbelow!(M)
                T = zeros(L, K)
                for k in 1:K
                    for l in 1:L
                        T[l, k] = rand(GeneralizedInverseGaussian(2ν, M[l, k], a[l] - 1))
                    end
                end
                ξ = T ./ sum(T, 2)
                ξ
            end
    )

    Gibbs_τ = Sampler([:τ],
        (τ, ξ, Λ, a, ν, L, K) ->
            begin
                S = 2vec(sum(abs(Λ) ./ ξ, 2))
                for l in 1:L
                    τ[l] = rand(GeneralizedInverseGaussian(2ν, S[l], K * a[l] - K))
                end
                τ
            end
    )

#   Gibbs_ξ = Sampler([:ξ],
#       (ξ, Λ, a, ν, L, K) ->
#           begin
#               A = a .- 1.0
#               @rput A
#               B = 2ν
#               M = 2abs(Λ)
#               boundbelow!(M)
#               @rput M
#               T = rcopy(reval("""
#                   Tmat <- matrix(0, nrow=$L, ncol=$K)
#                   for (l in 1:$L) {
#                       for(k in 1:$K) {
#                           Tmat[l, k] <- rgig(n=1, lambda=A[l], chi=M[l, k], psi=$B)
#                       }
#                   }
#                   Tmat
#                   """))
#               ξ = T ./ sum(T, 2)
#               ξ
#           end
#   )
#
#   Gibbs_τ = Sampler([:τ],
#       (ξ, Λ, a, ν, L, K) ->
#           begin
#               A = K .* a .- K
#               @rput A
#               B = 2ν
#               S = 2vec(sum(abs(Λ) ./ ξ, 2))
#               @rput S
#               τ = rcopy(reval("""
#                   tau <- vector('double', length=$L)
#                   for (l in 1:$L) {
#                       tau[l] <- rgig(n=1, lambda=A[l], chi=S[l], psi=$B)
#                   }
#                   tau
#                   """))
#               τ
#           end
#   )

    Gibbs_ν = Sampler([:ν],
        (τ, K, a, ν) ->
            begin
                u = K * sum(a) + shape(ν.distr)
                v = 1 / (sum(τ) + scale(ν.distr))
                rand(Gamma(u, v))
            end
    )

    Sample_a = Sampler([:idx_a],
        (idx_a, arr_a, ξ, τ, ν, L, K) ->
            begin
                dir = [Dirichlet(fill(a, K)) for a in arr_a]
                gam = [Gamma(K * a, 1 / ν) for a in arr_a]
                lp = zeros(arr_a)
                for l in 1:L
                    xi = ξ[l, :]
                    tau = τ[l]
                    for i in eachindex(lp, dir, gam)
                        lp[i] = logpdf(dir[i], xi) + logpdf(gam[i], tau)
                    end
                    lp .-= maximum(lp)
                    idx_a[l] = StatsBase.sample(1:length(arr_a), StatsBase.WeightVec(exp(lp)))
                end
                idx_a
            end
    )

    Gibbs_F = Sampler([:F],
        (θ, θ_var, μ, Λ, F, F_mean, b_var, L, N) ->
            begin
                A = Λ' ./ θ_var
                Sigma = inv(Λ * A + speye(L))
                mu = Sigma * ((θ .- μ') * A + F_mean)'
                for i in 1:N
                    F[i, :] = rand(MvNormal(mu[:, i], Hermitian(Sigma)))
                end
                F
            end
    )

    Gibbs_g = Sampler([:g],
        (F, Xb, g, g_var, Z, sumZg, q, L) ->
            begin
                A = F - Xb - sumZg
                for j in 1:size(Z, 2)
                    z = Z[:, j]
                    A += g[:, z]'
                    for r in unique(z)
                        in_block = z .== r
                        B = vec(sum(A[in_block, :], 1))
                        Sigma = 1 / (sum(in_block) + 1 / g_var[j])
                        mu = Sigma .* B
                        for l in 1:L
                            g[l, r] = rand(Normal(mu[l], sqrt(Sigma)))
                        end
                    end
                    A -= g[:, z]'
                end
                g
            end
    )

    Gibbs_b_full = Sampler([:b_full],
        (F, sumZg, b_full, b_var, X, ω, L, p) ->
            begin
                A = F - sumZg
                for l in 1:L
                    Xtilde = X * diagm(ω[l, :])
                    Sigma = inv(Xtilde' * Xtilde + diagm(1 ./ b_var))
                    mu = Sigma * (Xtilde' * A[:, l])
                    b_full[l, :] = rand(MvNormal(mu, Hermitian(Sigma)))
                end
                b_full
            end
    )

    Gibbs_ω = Sampler([:ω],
        (ω, F, sumZg, X, Xb, b_full, b, π, p, L) ->
            begin
                D = F - sumZg - Xb
                for l in 1:L
                    d = D[:, l]
                    for j in 1:p
                        x = X[:, j]
                        d += x .* b[l, j]
                        lprob_in = log(π[j]) - 0.5sumabs2(d - x .* b_full[l, j])
                        lprob_out = log(1 - π[j]) - 0.5sumabs2(d)
                        diff = min(lprob_in - lprob_out, 10)
                        prob_in = exp(diff) / (1 + exp(diff))
                        ω[l, j] = rand(Bernoulli(prob_in))
                        d -= x .* (ω[l, j] * b_full[l, j])
                    end
                end
                ω
            end
    )

    Gibbs_π = Sampler([:π],
        (π, ω, L, p) ->
            begin
                S = vec(sum(ω, 1))
                c0, d0 = params(π.distr)
                c = c0 .+ S
                d = d0 + L .- S
                for j in 1:p
                    π[j] = rand(Beta(c[j], d[j]))
                end
                π
            end
    )

    Gibbs_μ_var = Sampler([:μ_var],
        (μ, μ_var, K) ->
            begin
                u = 0.5K + shape(μ_var.distr)
                v = 0.5sumabs2(μ) + scale(μ_var.distr)
                rand(InverseGamma(u, v))
            end
    )

    Gibbs_g_var = Sampler([:g_var],
        (g, g_var, L, q, blocking_factor) ->
            begin
                num_blocking_factor = length(g_var)
                u = fill(shape(g_var.distr), num_blocking_factor)
                v = fill(scale(g_var.distr), num_blocking_factor)
                G = sumabs2(g, 1)
                for r in 1:q
                    u[blocking_factor[r]] += 0.5L
                    v[blocking_factor[r]] += 0.5G[r]
                end
                for i in 1:num_blocking_factor
                    g_var[i] = rand(InverseGamma(u[i], v[i]))
                end
                g_var
            end
    )

    Gibbs_b_var = Sampler([:b_var],
        (b, b_var, ω, p) ->
            begin
                u = 0.5vec(sum(ω, 1)) .+ shape(b_var.distr)
                v = 0.5vec(sumabs2(b, 1)) .+ scale(b_var.distr)
                for j in 1:p
                    b_var[j] = rand(InverseGamma(u[j], v[j]))
                end
                b_var
            end
    )

    Gibbs_θ_var = Sampler([:θ_var],
        (θ, θ_mean, θ_var, N, K) ->
            begin
                u = 0.5N + shape(θ_var.distr)
                v = 0.5vec(sumabs2(θ - θ_mean, 1)) .+ scale(θ_var.distr)
                for k in 1:K
                    θ_var[k] = rand(InverseGamma(u, v[k]))
                end
                θ_var
            end
    )

    scheme = [HMC_θ, Gibbs_F, Gibbs_Λ, Gibbs_ψ, Gibbs_τ, Gibbs_ξ,  Gibbs_ν,
              Sample_a, Gibbs_μ, Gibbs_b_full, Gibbs_ω, Gibbs_π, Gibbs_g,
              Gibbs_θ_var, Gibbs_μ_var, Gibbs_b_var, Gibbs_g_var]
    return scheme
end


