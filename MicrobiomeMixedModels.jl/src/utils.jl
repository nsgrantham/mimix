
function squeezesum(A, dim)
    dropdims(sum(A, dims=dim), dims=dim)
end

function boundbelow!(x, tol=1e-10)
    for i in eachindex(x)
        if x[i] < tol
            x[i] = tol
        end
    end
end

function proportionalize(Y::Matrix, tol=1e-10)
    ϕ = Y ./ sum(Y, dims=2)
    for i in eachindex(ϕ)
        if ϕ[i] < tol
            ϕ[i] = tol
        end
    end
    ϕ ./ sum(ϕ, dims=2)
end

function clr(θ)
    x = exp.(θ)
    S = sum(x, dims=2)
    x ./ S
end

function iclr(ϕ)
    x = log.(ϕ)
    x .- mean(x, dims=2)
end

function logfgrad(θ, θ_mean, θ_var, y, m)
    e = exp.(θ)
    S = sum(e)
    resid = θ - θ_mean
    logf = dot(y, θ) - m * log(S) - 0.5sum((resid.^2) ./ θ_var)
    grad = y - m .* e ./ S - resid ./ θ_var
    logf, grad
end

function get_post(sim, data, param)
    N = data[:N]
    K = data[:K]
    L = data[:L]
    p = data[:p]
    num_blocking_factors = data[:num_blocking_factors]
    param_names = Dict{Symbol, Any}(
        :Y => ["Y[$i, $k]" for i in 1:N, k in 1:K],
        :ϕ => ["ϕ[$i, $k]" for i in 1:N, k in 1:K],
        :θ => ["θ[$i, $k]" for i in 1:N, k in 1:K],
        :θ_mean => ["θ_mean[$i, $k]" for i in 1:N, k in 1:K],
        :θ_var => ["θ_var[$k]" for k in 1:K],
        :μ => ["μ[$k]" for k in 1:K],
        :μ_var => ["μ_var"],
        :β => ["β[$k, $j]" for k in 1:K, j in 1:p],
        :β_full => ["β[$k, $j]" for k in 1:K, j in 1:p],
        :β_var => ["β_var[$j]" for j in 1:p],
        :γ => ["γ[$k, $r]" for k in 1:K, r in 1:num_blocking_factors],
        :γ_var => ["γ_var[$r]" for r in 1:num_blocking_factors],
        :ω => ["ω[$l, $j]" for j in 1:p, l in 1:L],
        :π => ["π[$j]" for j in 1:p],
        :Λ => ["Λ[$l, $k]" for l in 1:L, k in 1:K],
        :λ_var => ["λ_var[$l, $k]" for l in 1:L, k in 1:K],
        :a => ["a[$l]" for l in 1:L],
        :ψ => ["ψ[$l, $k]" for l in 1:L, k in 1:K],
        :ξ => ["ξ[$l, $k]" for l in 1:L, k in 1:K],
        :τ => ["τ[$l]" for l in 1:L],
        :ν => ["ν"],
        :g => ["g[$l, $r]" for l in 1:L, r in 1:num_blocking_factors],
        :g_var => ["g_var[$r]" for r in 1:num_blocking_factors],
        :b => ["b[$l, $j]" for l in 1:L, j in 1:p],
        :b_full => ["b[$l, $j]" for l in 1:L, j in 1:p],
        :b_var => ["b_var[$j]" for j in 1:p],
    )
    _, _, chains = size(sim)
    names = param_names[param]
    post = vcat([sim[:, vec(names), c].value for c in 1:chains]...)
    post = DataFrame(dropdims(post, dims=3))
    names!(post, [Symbol(name) for name in vec(names)])
    return post
end

function parse_config(conf)
    parsed_conf = Dict{Symbol, Any}()
    for key in keys(conf)
        param = key
        for (name, letter) in greek
            if startswith(param, name)
                param = replace(param, name => letter)
            end
        end
        parsed_conf[Symbol(param)] = conf[key]
    end
    return parsed_conf
end
    
greek = Dict(
    "beta" => "β",
    "epsilon" => "ϵ",
    "gamma" => "γ",
    "lambda" => "λ",
    "Lambda" => "Λ",
    "mu" => "μ",
    "nu" => "ν",
    "omega" => "ω",
    "phi" => "ϕ",
    "pi" => "π",
    "psi" => "ψ",  # must come after epsilon, or will change epsilon -> eψlon :(
    "tau" => "τ",
    "theta" => "θ",
    "xi" => "ξ"
)

latin = Dict([(v, k) for (k, v) in greek])

function load_config(filename::AbstractString)
    @assert isfile(filename)
    conf = YAML.load(open(filename, "r"))
    return parse_config(conf)
end

load_config(filenames::Vector{T}) where T <: AbstractString = merge!([load_config(f) for f in filenames]...)

