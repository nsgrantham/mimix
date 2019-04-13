
function get_inits(::MIMIXGaussian, params::Dict{Symbol, Any}, data::Dict{Symbol, Any})
    θ_init = iclr(proportionalize(data[:Y]))
    ω_init = params[:ω] == "ones" ? ones(data[:L], data[:p]) : zeros(data[:L], data[:p])
    π_init = params[:π] == "even" ? [0.5 for _ in 1:data[:p]] : params[:π]
    b_init = params[:b] == "zeros" ? zeros(data[:L], data[:p]) : params[:b]
    g_init = params[:g] == "zeros" ? zeros(data[:L], data[:q]) : params[:g]
    b_var_init = params[:b_var] == "ones" ? ones(data[:p]) : params[:b_var]
    g_var_init = params[:g_var] == "ones" ? ones(data[:num_blocking_factors]) : params[:g_var]
    
    inits = Dict{Symbol, Any}(
        :Y => data[:Y],
        :θ => θ_init,
        :θ_var => vec(var(θ_init, dims=1)),
        :μ => vec(mean(θ_init, dims=1)),
        :μ_var => params[:μ_var],
        :Λ => Matrix(1.0I, data[:L], data[:K]),
        :λ_var => ones(data[:L]),
        :F => zeros(data[:N], data[:L]),
        :b_full => b_init,
        :g => g_init,
        :ω => ω_init,
        :π => π_init,
        :b_var => b_var_init,
        :g_var => g_var_init
    )

    return inits
end

function get_model(::MIMIXGaussian, monitor::Dict{Symbol, Any}, hyper::Dict{Symbol, Any})
    
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
            (μ, ΛF) -> ΛF .+ transpose(μ),
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

        # Full-dimensional fixed effect estimates (monitor these for OTU-level inference!)

        β = Logical(2,
            (Λ, b) -> transpose(Λ) * b,
            monitor[:β]
        ),

        # Factor analysis w/ loadings matrix Λ (L by K) and scores F (N by L)

        ΛF = Logical(2,
            (Λ, F) -> F * Λ,
            false
        ),

        Λ = Stochastic(2,
            (λ_var, L, K) -> MultivariateDistribution[
                MvNormal(zeros(K), sqrt(λ_var[l])) for l in 1:L
            ],
            monitor[:Λ]
        ),

        λ_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            monitor[:λ_var]
        ),

        F = Stochastic(2,
            (F_mean, N) -> MultivariateDistribution[
                MvNormal(F_mean[i, :], 1.0) for i in 1:N
            ],
            false
        ),

        F_mean = Logical(2,
            (Xb, sumZg) -> Xb + sumZg,
            false
        ),

        # Blocking factor random effects

        sumZg = Logical(2,
            (Zg) -> transpose(squeezesum(Zg, 3)),
            false
        ),

        Zg = Logical(3,
            (Z, g) -> g[:, Z],
            false
        ),

        g = Stochastic(2,
            (q, L, g_var, blocking_factor) -> UnivariateDistribution[
                Normal(0.0, sqrt(g_var[blocking_factor[r]])) for l in 1:L, r in 1:q
            ],
            monitor[:g]
        ),

        g_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            monitor[:g_var]
        ),

        # Fixed treatment effects w/ spike-and-slab variable selection

        Xb = Logical(2,
            (X, b) -> X * transpose(b),
            false
        ),

        b = Logical(2,
            (b_full, ω) -> ω .* b_full,
            monitor[:b]
        ),

        b_full = Stochastic(2,
            (p, L, b_var) -> MultivariateDistribution[
                MvNormal(zeros(p), sqrt.(b_var)) for l in 1:L
            ],
            monitor[:b_full]
        ),

        b_var = Stochastic(1,
            () -> InverseGamma(1, 1),
            monitor[:b_var]
        ),

        ω = Stochastic(2,
            (p, L, π) -> UnivariateDistribution[
                Bernoulli(π[j]) for l in 1:L, j in 1:p
            ],
            monitor[:ω]
        ),

        π = Stochastic(1,
            (L) -> Beta(1, L),
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
                    v1 += 0.5eps * grad0
                    for step in 1:hyper[:hmc_steps]
                        u1 += eps * v1
                        logf1, grad1 = logfgrad(u1, θi_mean, θ_var, yi, mi)
                        v1 += eps * grad1
                    end
                    v1 -= 0.5eps * grad1
                    v1 *= -1.0
                    Kv0 = 0.5sum(abs2, v0)
                    Kv1 = 0.5sum(abs2, v1)
                    if rand() < exp((logf1 - Kv1) - (logf0 - Kv0))
                        acc += 1
                        θ[i, :] = u1
                    end
                end
                if hyper[:hmc_verbose]
                    println("HMC accept rate: $(round(100acc/N, 1))%, ")
                end
                θ
            end
        ),

        Sampler(:μ, (θ, θ_var, μ, μ_var, ΛF, N, K) ->
            begin
                Sigma = 1 ./ (N ./ θ_var .+ 1 / μ_var)
                mu = Sigma .* vec(sum(θ - ΛF, dims=1)) ./ θ_var
                for k in 1:K
                    μ[k] = rand(Normal(mu[k], sqrt(Sigma[k])))
                end
                μ
            end
        ),

        Sampler(:Λ, (θ, θ_var, μ, Λ, λ_var, F, K, L) ->
            begin
                FtF = transpose(F) * F
                B = (transpose(F) * (θ .- transpose(μ))) ./ transpose(θ_var)
                Sigma = zeros(L, L)
                for k in 1:K
                    Sigma[:, :] = inv(cholesky(Hermitian(Diagonal(1 ./ λ_var) + FtF ./ θ_var[k])))
                    Λ[:, k] = rand(MvNormal(Sigma * B[:, k], Sigma))
                end
                Λ
            end
        ),

        Sampler(:F, (θ, θ_var, μ, Λ, F, F_mean, b_var, L, N) ->
            begin
                A = transpose(Λ) ./ θ_var
                Sigma = inv(cholesky(Hermitian(Λ * A + Matrix(1.0I, L, L))))
                mu = Sigma * transpose((θ .- transpose(μ)) * A + F_mean)
                for i in 1:N
                    F[i, :] = rand(MvNormal(mu[:, i], Sigma))
                end
                F
            end
        ),

        Sampler(:g, (F, Xb, g, g_var, Z, sumZg, q, L) ->
            begin
                A = F - Xb - sumZg
                for j in 1:size(Z, 2)
                    z = Z[:, j]
                    A += transpose(g[:, z])
                    for r in unique(z)
                        in_block = z .== r
                        B = vec(sum(A[in_block, :], dims=1))
                        Sigma = 1 / (sum(in_block) + 1 / g_var[j])
                        mu = Sigma .* B
                        for l in 1:L
                            g[l, r] = rand(Normal(mu[l], sqrt(Sigma)))
                        end
                    end
                    A -= transpose(g[:, z])
                end
                g
            end
        ),

        Sampler(:b_full, (F, sumZg, b_full, b_var, X, ω, L, p) ->
            begin
                A = F - sumZg
                for l in 1:L
                    Xtilde = X * Diagonal(ω[l, :])
                    Sigma = inv(cholesky(Hermitian(transpose(Xtilde) * Xtilde + Diagonal(1 ./ b_var))))
                    mu = Sigma * (transpose(Xtilde) * A[:, l])
                    b_full[l, :] = rand(MvNormal(mu, Sigma))
                end
                b_full
            end
        ),

        Sampler(:ω, (ω, F, sumZg, X, Xb, b_full, b, π, p, L) ->
            begin
                D = F - sumZg - Xb
                for l in 1:L
                    d = D[:, l]
                    for j in 1:p
                        x = X[:, j]
                        d += x .* b[l, j]
                        lprob_in = log(π[j]) - 0.5sum(abs2, d - x .* b_full[l, j])
                        lprob_out = log(1 - π[j]) - 0.5sum(abs2, d)
                        diff = min(lprob_in - lprob_out, 10)
                        prob_in = exp(diff) / (1 + exp(diff))
                        ω[l, j] = rand(Bernoulli(prob_in))
                        d -= x .* (ω[l, j] * b_full[l, j])
                    end
                end
                ω
            end
        ),

        Sampler(:π, (π, ω, L, p) ->
            begin
                S = vec(sum(ω, dims=1))
                c0, d0 = params(π.distr)
                c = c0 .+ S
                d = d0 + L .- S
                for j in 1:p
                    π[j] = rand(Beta(c[j], d[j]))
                end
                π
            end
        ),

        Sampler(:μ_var, (μ, μ_var, K) ->
            begin
                u = 0.5K + shape(μ_var.distr)
                v = 0.5sum(abs2, μ) + scale(μ_var.distr)
                rand(InverseGamma(u, v))
            end
        ),

        Sampler(:g_var, (g, g_var, L, q, blocking_factor) ->
            begin
                num_blocking_factor = length(g_var)
                u = fill(shape(g_var.distr), num_blocking_factor)
                v = fill(scale(g_var.distr), num_blocking_factor)
                G = sum(abs2, g, dims=1)
                for r in 1:q
                    u[blocking_factor[r]] += 0.5L
                    v[blocking_factor[r]] += 0.5G[r]
                end
                for i in 1:num_blocking_factor
                    g_var[i] = rand(InverseGamma(u[i], v[i]))
                end
                g_var
            end
        ),

        Sampler(:b_var, (b, b_var, ω, p) ->
            begin
                u = 0.5vec(sum(ω, dims=1)) .+ shape(b_var.distr)
                v = 0.5vec(sum(abs2, b, dims=1)) .+ scale(b_var.distr)
                for j in 1:p
                    b_var[j] = rand(InverseGamma(u[j], v[j]))
                end
                b_var
            end
        ),

        Sampler(:θ_var, (θ, θ_mean, θ_var, N, K) ->
            begin
                u = 0.5N + shape(θ_var.distr)
                v = 0.5vec(sum(abs2, θ - θ_mean, dims=1)) .+ scale(θ_var.distr)
                for k in 1:K
                    θ_var[k] = rand(InverseGamma(u, v[k]))
                end
                θ_var
            end
        ),

        Sampler(:λ_var, (λ_var, Λ, K, L) ->
        begin
            u = 0.5K + shape(λ_var.distr)
            v = 0.5vec(sum(abs2, Λ, dims=2)) .+ scale(λ_var.distr)
            [rand(InverseGamma(u, v[l])) for l in 1:L]
        end
    )
    ]

    setsamplers!(model, samplers)
    return model
end


