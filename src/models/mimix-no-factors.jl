
immutable MIMIXNoFactors end

function fit(mm::MIMIXNoFactors, Y, X, Z;
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

function datadict(mm::MIMIXNoFactors, Y, X, Z)
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
    return d
end

function inits(::MIMIXNoFactors, d::Dict{Symbol, Any})
    θ_init = iclr(proportionalize(d[:Y]))
    i = [
        Dict{Symbol, Any}(
            :Y => d[:Y],
            :θ => θ_init,
            :θ_var => vec(var(θ_init, 1)),
            :μ => vec(mean(θ_init, 1)),
            :μ_var => 1.0,
            :ω => ones(d[:K], d[:p]),
            :π => [0.5 for j in 1:d[:p]],
            :β_full => zeros(d[:K], d[:p]),
            :γ => zeros(d[:K], d[:q]),
            :β_var => ones(d[:p]),
            :γ_var => ones(d[:num_blocking_factor])
        )
    ]
    return i
end

function model(::MIMIXNoFactors, monitor)
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
            (μ, Xβ, sumZγ) -> Xβ + sumZγ .+ μ',
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

        # Blocking factor random effects

        sumZγ = Logical(2,
            (Zγ) -> squeezesum(Zγ, 3)',
            :sumZγ in monitor
        ),

        Zγ = Logical(3,
            (Z, γ) -> γ[:, Z],
            :Zγ in monitor
        ),

        γ = Stochastic(2,
            (q, K, γ_var, blocking_factor) -> UnivariateDistribution[
                Normal(0.0, sqrt(γ_var[blocking_factor[r]])) for k in 1:K, r in 1:q
            ],
            :γ in monitor
        ),

        γ_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            :γ_var in monitor
        ),

        # Fixed treatment effects w/ spike-and-slab variable selection

        Xβ = Logical(2,
            (X, β) -> X * β',
            :Xβ in monitor
        ),

        β = Logical(2,
            (β_full, ω) -> ω .* β_full,
            :β in monitor
        ),

        β_full = Stochastic(2,
            (p, K, β_var) -> MultivariateDistribution[
                MvNormal(zeros(p), sqrt(β_var)) for k in 1:K
            ],
            :β_full in monitor
        ),

        β_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            :β_var in monitor
        ),

        ω = Stochastic(2,
            (p, K, π) -> UnivariateDistribution[
                Bernoulli(π[j]) for k in 1:K, j in 1:p
            ],
            :ω in monitor
        ),

        π = Stochastic(1,
            (K) -> Beta(1, K),
            :π in monitor
        )
    )
end

function samplers(::MIMIXNoFactors, epsilon::Float64, steps::Int, hmc_verbose::Bool)

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
                    v1 += 0.5 * eps * grad0
                    for step in 1:steps
                        u1 += eps * v1
                        logf1, grad1 = logfgrad(u1, θi_mean, θ_var, yi, mi)
                        v1 += eps * grad1
                    end
                    v1 -= 0.5 * eps * grad1
                    v1 *= -1.0
                    Kv0 = 0.5 * sumabs2(v0)
                    Kv1 = 0.5 * sumabs2(v1)
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
        (θ, θ_var, μ, μ_var, Xβ, sumZγ, N, K) ->
            begin
                Sigma = 1 ./ (N ./ θ_var .+ 1 / μ_var)
                mu = Sigma .* vec(sum(θ - Xβ - sumZγ, 1)) ./ θ_var
                for k in 1:K
                    μ[k] = rand(Normal(mu[k], sqrt(Sigma[k])))
                end
                μ
            end
    )

    Gibbs_γ = Sampler([:γ],
        (θ, θ_var, Xβ, μ, γ, γ_var, Z, sumZγ, q, K) ->
            begin
                A = θ - Xβ - sumZγ .- μ'
                for j in 1:size(Z, 2)
                    z = Z[:, j]
                    A += γ[:, z]'
                    for r in unique(z)
                        in_block = z .== r
                        B = vec(sum(A[in_block, :], 1))
                        Sigma = 1 ./ (sum(in_block) ./ θ_var .+ 1 / γ_var[j])
                        mu = Sigma .* B ./ θ_var
                        for k in 1:K
                            γ[k, r] = rand(Normal(mu[k], sqrt(Sigma[k])))
                        end
                    end
                    A -= γ[:, z]'
                end
                γ
            end
    )

    Gibbs_β_full = Sampler([:β_full],
        (θ, θ_var, sumZγ, μ, β_full, β_var, X, ω, K, p) ->
            begin
                A = θ - sumZγ .- μ'
                for k in 1:K
                    Xtilde = X * diagm(ω[k, :])
                    Sigma = inv(Xtilde' * Xtilde ./ θ_var[k] + diagm(1 ./ β_var))
                    mu = Sigma * ((Xtilde' * A[:, k]) ./ θ_var[k])
                    β_full[k, :] = rand(MvNormal(mu, Hermitian(Sigma)))
                end
                β_full
            end
    )

    Gibbs_ω = Sampler([:ω],
        (ω, θ, θ_var, sumZγ, X, Xβ, μ, β, β_full, π, p, K) ->
            begin
                D = θ - sumZγ - Xβ .- μ'
                for k in 1:K
                    d = D[:, k]
                    for j in 1:p
                        x = X[:, j]
                        d += x .* β[k, j]
                        lprob_in = log(π[j]) - 0.5sumabs2(d - x .* β_full[k, j]) / θ_var[k]
                        lprob_out = log(1 - π[j]) - 0.5sumabs2(d) / θ_var[k]
                        diff = min(lprob_in - lprob_out, 10.0)
                        prob_in = exp(diff) / (1 + exp(diff))
                        ω[k, j] = rand(Bernoulli(prob_in))
                        d -= x .* (ω[k, j] * β_full[k, j])
                    end
                end
                ω
            end
    )

    Gibbs_π = Sampler([:π],
        (π, ω, K, p) ->
            begin
                S = vec(sum(ω, 1))
                c0, d0 = params(π.distr)
                c = c0 .+ S
                d = d0 + K .- S
                [rand(Beta(c[j], d[j])) for j in 1:p]
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

    Gibbs_γ_var = Sampler([:γ_var],
        (γ, γ_var, K, q, blocking_factor) ->
            begin
                num_blocking_factor = length(γ_var)
                u = fill(shape(γ_var.distr), num_blocking_factor)
                v = fill(scale(γ_var.distr), num_blocking_factor)
                G = sumabs2(γ, 1)
                for r in 1:q
                    u[blocking_factor[r]] += 0.5K
                    v[blocking_factor[r]] += 0.5G[r]
                end
                for i in 1:num_blocking_factor
                    γ_var[i] = rand(InverseGamma(u[i], v[i]))
                end
                γ_var
            end
    )

    Gibbs_β_var = Sampler([:β_var],
        (β, β_var, ω, p) ->
            begin
                u = 0.5vec(sum(ω, 1)) .+ shape(β_var.distr)
                v = 0.5vec(sumabs2(β, 1)) .+ scale(β_var.distr)
                [rand(InverseGamma(u[j], v[j])) for j in 1:p]
            end
    )

    Gibbs_θ_var = Sampler([:θ_var],
        (θ, θ_mean, θ_var, N, K) ->
            begin
                u = 0.5N + shape(θ_var.distr)
                v = 0.5vec(sumabs2(θ - θ_mean, 1)) .+ scale(θ_var.distr)
                [rand(InverseGamma(u, v[k])) for k in 1:K]
            end
    )

    scheme = [HMC_θ, Gibbs_μ, Gibbs_β_full, Gibbs_ω, Gibbs_π, Gibbs_γ,
              Gibbs_θ_var, Gibbs_μ_var, Gibbs_β_var, Gibbs_γ_var]
    return scheme
end
