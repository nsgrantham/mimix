
function get_inits(::MIMIXNoFactors, params::Dict{Symbol, Any}, data::Dict{Symbol, Any})
    θ_init = iclr(proportionalize(data[:Y]))
    ω_init = params[:ω] == "ones" ? ones(data[:K], data[:p]) : zeros(data[:K], data[:p])
    π_init = params[:π] == "even" ? [0.5 for _ in 1:data[:p]] : params[:π]
    β_init = params[:β] == "zeros" ? zeros(data[:K], data[:p]) : params[:β]
    γ_init = params[:γ] == "zeros" ? zeros(data[:K], data[:q]) : params[:γ]
    β_var_init = params[:β_var] == "ones" ? ones(data[:p]) : params[:β_var]
    γ_var_init = params[:γ_var] == "ones" ? ones(data[:num_blocking_factors]) : params[:γ_var]
    
    inits = Dict{Symbol, Any}(
        :Y => data[:Y],
        :θ => θ_init,
        :θ_var => vec(var(θ_init, dims=1)),
        :μ => vec(mean(θ_init, dims=1)),
        :μ_var => params[:μ_var],
        :ω => ω_init,
        :π => π_init,
        :β_full => β_init,
        :γ => γ_init,
        :β_var => β_var_init,
        :γ_var => γ_var_init
    )
    return inits
end

function get_model(::MIMIXNoFactors, monitor::Dict{Symbol, Any}, hyper::Dict{Symbol, Any})
    
    model = Model(
        # High-dimensional counts
        Y = Stochastic(2,
            (m, ϕ, N) -> MultivariateDistribution[
                Multinomial(m[i], ϕ[i, :]) for i in 1:N
            ],
            monitor[:Y]
        ),

        ϕ = Logical(2,
            (θ) -> clr(θ),
            monitor[:ϕ]
        ),

        # Sample-specific random effects

        θ = Stochastic(2,
            (θ_mean, θ_var, N) -> MultivariateDistribution[
                MvNormal(θ_mean[i, :], sqrt.(θ_var)) for i in 1:N
            ],
            monitor[:θ]
        ),

        θ_mean = Logical(2,
            (μ, Xβ, sumZγ) -> Xβ + sumZγ .+ transpose(μ),
            monitor[:θ_mean]
        ),

        θ_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            monitor[:θ_var]
        ),

        # Population mean

        μ = Stochastic(1,
            (K, μ_var) -> MvNormal(K, sqrt(μ_var)),
            monitor[:μ]
        ),

        μ_var = Stochastic(
            () -> InverseGamma(1, 1),
            monitor[:μ_var]
        ),

        # Blocking factor random effects

        sumZγ = Logical(2,
            (Zγ) -> transpose(squeezesum(Zγ, 3)),
            false
        ),

        Zγ = Logical(3,
            (Z, γ) -> γ[:, Z],
            false
        ),

        γ = Stochastic(2,
            (q, K, γ_var, blocking_factor) -> UnivariateDistribution[
                Normal(0.0, sqrt(γ_var[blocking_factor[r]])) for k in 1:K, r in 1:q
            ],
            monitor[:γ]
        ),

        γ_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            monitor[:γ_var]
        ),

        # Fixed treatment effects w/ spike-and-slab variable selection

        Xβ = Logical(2,
            (X, β) -> X * transpose(β),
            false
        ),

        β = Logical(2,
            (β_full, ω) -> ω .* β_full,
            monitor[:β]
        ),

        β_full = Stochastic(2,
            (p, K, β_var) -> MultivariateDistribution[
                MvNormal(zeros(p), sqrt.(β_var)) for k in 1:K
            ],
            monitor[:β_full]
        ),

        β_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            monitor[:β_var]
        ),

        ω = Stochastic(2,
            (p, K, π) -> UnivariateDistribution[
                Bernoulli(π[j]) for k in 1:K, j in 1:p
            ],
            monitor[:ω]
        ),

        π = Stochastic(1,
            (K) -> Beta(1, K),
            monitor[:π]
        )
    )

    samplers = [
        Sampler(:θ, (Y, m, θ, θ_mean, θ_var, N, K) ->
            begin
                acc = 0
                for i in 1:N
                    eps = rand(Exponential(hyper[:hmc_epsilon]))
                    yi = Y[i, :]
                    mi = m[i]
                    u1 = θ[i, :]
                    θi_mean = θ_mean[i, :]
                    logf0, grad0 = logf1, grad1 = logfgrad(u1, θi_mean, θ_var, yi, mi)
                    v0 = v1 = randn(length(u1))
                    v1 += 0.5 * eps * grad0
                    for step in 1:hyper[:hmc_steps]
                        u1 += eps * v1
                        logf1, grad1 = logfgrad(u1, θi_mean, θ_var, yi, mi)
                        v1 += eps * grad1
                    end
                    v1 -= 0.5 * eps * grad1
                    v1 *= -1.0
                    Kv0 = 0.5 * sum(abs2, v0)
                    Kv1 = 0.5 * sum(abs2, v1)
                    if rand() < exp((logf1 - Kv1) - (logf0 - Kv0))
                        acc += 1
                        θ[i, :] = u1
                    end
                end
                if hyper[:hmc_verbose]
                    println("HMC accept rate: $(round(100acc/N, digits=1))%, ")
                end
                θ
            end
        ),

        Sampler(:μ, (θ, θ_var, μ, μ_var, Xβ, sumZγ, N, K) ->
            begin
                Sigma = 1 ./ (N ./ θ_var .+ 1 / μ_var)
                mu = Sigma .* vec(sum(θ - Xβ - sumZγ, dims=1)) ./ θ_var
                for k in 1:K
                    μ[k] = rand(Normal(mu[k], sqrt(Sigma[k])))
                end
                μ
            end
        ),

        Sampler([:γ], (θ, θ_var, Xβ, μ, γ, γ_var, Z, sumZγ, q, K) ->
            begin
                A = θ - Xβ - sumZγ .- transpose(μ)
                for j in 1:size(Z, 2)
                    z = Z[:, j]
                    A += transpose(γ[:, z])
                    for r in unique(z)
                        in_block = z .== r
                        B = vec(sum(A[in_block, :], dims=1))
                        Sigma = 1 ./ (sum(in_block) ./ θ_var .+ 1 / γ_var[j])
                        mu = Sigma .* B ./ θ_var
                        for k in 1:K
                            γ[k, r] = rand(Normal(mu[k], sqrt(Sigma[k])))
                        end
                    end
                    A -= transpose(γ[:, z])
                end
                γ
            end
        ),

        Sampler(:β_full, (θ, θ_var, sumZγ, μ, β_full, β_var, X, ω, K, p) ->
            begin
                A = θ - sumZγ .- transpose(μ)
                for k in 1:K
                    Xtilde = X * Diagonal(ω[k, :])
                    Sigma = inv(cholesky(Hermitian(transpose(Xtilde) * Xtilde ./ θ_var[k] + Diagonal(1 ./ β_var))))
                    mu = Sigma * ((transpose(Xtilde) * A[:, k]) ./ θ_var[k])
                    β_full[k, :] = rand(MvNormal(mu, Sigma))
                end
                β_full
            end
        ),

        Sampler(:ω, (ω, θ, θ_var, sumZγ, X, Xβ, μ, β, β_full, π, p, K) ->
            begin
                D = θ - sumZγ - Xβ .- transpose(μ)
                for k in 1:K
                    d = D[:, k]
                    for j in 1:p
                        x = X[:, j]
                        d += x .* β[k, j]
                        lprob_in = log(π[j]) - 0.5sum(abs2, d - x .* β_full[k, j]) / θ_var[k]
                        lprob_out = log(1 - π[j]) - 0.5sum(abs2, d) / θ_var[k]
                        diff = min(lprob_in - lprob_out, 10.0)
                        prob_in = exp(diff) / (1 + exp(diff))
                        ω[k, j] = rand(Bernoulli(prob_in))
                        d -= x .* (ω[k, j] * β_full[k, j])
                    end
                end
                ω
            end
        ),

        Sampler(:π, (π, ω, K, p) ->
            begin
                S = vec(sum(ω, dims=1))
                c0, d0 = params(π.distr)
                c = c0 .+ S
                d = d0 + K .- S
                [rand(Beta(c[j], d[j])) for j in 1:p]
            end
        ),

        Sampler(:μ_var, (μ, μ_var, K) ->
            begin
                u = 0.5K + shape(μ_var.distr)
                v = 0.5sum(abs2, μ) + scale(μ_var.distr)
                rand(InverseGamma(u, v))
            end
        ),

        Sampler(:γ_var, (γ, γ_var, K, q, blocking_factor) ->
            begin
                num_blocking_factor = length(γ_var)
                u = fill(shape(γ_var.distr), num_blocking_factor)
                v = fill(scale(γ_var.distr), num_blocking_factor)
                G = sum(abs2, γ, dims=1)
                for r in 1:q
                    u[blocking_factor[r]] += 0.5K
                    v[blocking_factor[r]] += 0.5G[r]
                end
                for i in 1:num_blocking_factor
                    γ_var[i] = rand(InverseGamma(u[i], v[i]))
                end
                γ_var
            end
        ),

        Sampler(:β_var, (β, β_var, ω, p) ->
            begin
                u = 0.5vec(sum(ω, dims=1)) .+ shape(β_var.distr)
                v = 0.5vec(sum(abs2, β, dims=1)) .+ scale(β_var.distr)
                [rand(InverseGamma(u[j], v[j])) for j in 1:p]
            end
        ),

        Sampler(:θ_var, (θ, θ_mean, θ_var, N, K) ->
            begin
                u = 0.5N + shape(θ_var.distr)
                v = 0.5vec(sum(abs2, θ - θ_mean, dims=1)) .+ scale(θ_var.distr)
                [rand(InverseGamma(u, v[k])) for k in 1:K]
            end
        )
    ]

    setsamplers!(model, samplers)
    return model
end
