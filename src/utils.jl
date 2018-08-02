
function squeezesum(A, dim)
    squeeze(sum(A, dim), dim)
end

function boundbelow!(x, tol=1e-10)
    for i in eachindex(x)
        if x[i] < tol
            x[i] = tol
        end
    end
end

function proportionalize(Y::Matrix, tol=1e-10)
    ϕ = Y ./ sum(Y, 2)
    for i in eachindex(ϕ)
        if ϕ[i] < tol
            ϕ[i] = tol
        end
    end
    ϕ ./ sum(ϕ, 2)
end

function clr(θ)
    x = exp.(θ)
    S = sum(x, 2)
    x ./ S
end

function iclr(ϕ)
    x = log.(ϕ)
    x .- mean(x, 2)
end

function logfgrad(θ, θ_mean, θ_var, y, m)
    e = exp.(θ)
    S = sum(e)
    resid = θ - θ_mean
    logf = dot(y, θ) - m * log(S) - 0.5sum((resid.^2) ./ θ_var)
    grad = y - m .* e ./ S - resid ./ θ_var
    logf, grad
end

function parse_config(conf)
    greek = Dict(
        "beta" => "β",
        "epsilon" => "ϵ",
        "gamma" => "γ",
        "lambda" => "λ",
        "mu" => "μ",
        "nu" => "ν",
        "omega" => "ω",
        "phi" => "ϕ",
        "pi" => "π",
        "psi" => "ψ",  # must come after epsilon
        "tau" => "τ",
        "theta" => "θ",
        "xi" => "ξ"
    )
    parsed_conf = Dict{Symbol, Any}()
    for key in keys(conf)
        param = key
        for (name, letter) in greek
            if startswith(param, name)
                param = replace(param, name, letter)
            end
        end
        parsed_conf[Symbol(param)] = conf[key]
    end
    return parsed_conf
end

function load_config(filename::AbstractString)
    @assert isfile(filename)
    conf = YAML.load(open(filename, "r"))
    return parse_config(conf)
end

load_config(filenames::Vector{T}) where T <: AbstractString = merge!([load_config(f) for f in filenames]...)