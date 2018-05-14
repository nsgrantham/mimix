
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
    x = exp(θ)
    S = sum(x, 2)
    x ./ S
end

function iclr(ϕ)
    x = log(ϕ)
    x .- mean(x, 2)
end

function logfgrad(θ, θ_mean, θ_var, y, m)
    expθ = exp(θ)
    S = sum(expθ)
    resid = θ - θ_mean
    logf = dot(y, θ) - m * log(S) - 0.5sum((resid.^2) ./ θ_var)
    grad = y - m .* expθ ./ S - resid ./ θ_var
    logf, grad
end
