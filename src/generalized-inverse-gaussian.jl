
struct GeneralizedInverseGaussian{T<:Real} <: ContinuousUnivariateDistribution
    a::T
    b::T
    p::T

    function GeneralizedInverseGaussian{T}(a::T, b::T, p::T) where {T}
        @assert a > zero(a) && b > zero(b)
        new{T}(a, b, p)
    end
end

GeneralizedInverseGaussian(a::T, b::T, p::T) where {T<:Real} = GeneralizedInverseGaussian{T}(a, b, p)
GeneralizedInverseGaussian(a::Real, b::Real, p::Real) = GeneralizedInverseGaussian(promote(a, b, p)...)
GeneralizedInverseGaussian(a::Integer, b::Integer, p::Integer) = GeneralizedInverseGaussian(Float64(a), Float64(b), Float64(p))

Distributions.@distr_support GeneralizedInverseGaussian 0.0 Inf


#### Conversions

function convert(::Type{GeneralizedInverseGaussian{T}}, a::S, b::S, p::S) where {T <: Real, S <: Real}
    GeneralizedInverseGaussian(T(a), T(b), T(p))
end
function convert(::Type{GeneralizedInverseGaussian{T}}, d::GeneralizedInverseGaussian{S}) where {T <: Real, S <: Real}
    GeneralizedInverseGaussian(T(d.a), T(d.b), T(d.p))
end

#### Parameters

params(d::GeneralizedInverseGaussian) = (d.a, d.b, d.p)
@inline partype(d::GeneralizedInverseGaussian{T}) where {T <: Real} = T


#### Statistics

function mean(d::GeneralizedInverseGaussian)
    a, b, p = params(d)
    q = sqrt(a * b)
    (sqrt(b) * besselk(p + 1, q)) / (sqrt(a) * besselk(p, q))
end

function var(d::GeneralizedInverseGaussian)
    a, b, p = params(d)
    q = sqrt(a * b)
    r = besselk(p, q)
    (b / a) * ((besselk(p + 2, q) / r) - (besselk(p + 1, q) / r)^2)
end

mode(d::GeneralizedInverseGaussian) = ((d.p - 1) + sqrt((d.p - 1)^2 + d.a * d.b)) / d.a


#### Evaluation

function pdf(d::GeneralizedInverseGaussian{T}, x::Real) where {T <: Real}
    if x > 0
        a, b, p = params(d)
        (((a / b)^(p / 2)) / (2 * besselk(p, sqrt(a * b)))) * (x^(p - 1)) * exp(- (a * x + b / x) / 2)
    else
        zero(T)
    end
end

function logpdf(d::GeneralizedInverseGaussian{T}, x::Real) where {T <: Real}
    if x > 0
        a, b, p = params(d)
        (p / 2) * (log(a) - log(b)) - log(2 * besselk(p, sqrt(a * b))) + (p - 1) * log(x) - (a * x + b / x) / 2
    else
        -T(Inf)
    end
end


function cdf(d::GeneralizedInverseGaussian{T}, x::Real) where {T <: Real}
    if x > 0
        # See eq. (5) in Lemonte & Cordeiro (2011) 
        # Statistics & Probability Letters 81:506–517
        # F(x) = 1 - (ρ + σ), where ρ and σ are infinite sums
        # calculated up to truncation below
        a, b, p = params(d)
        c = (((a / b)^(p / 2)) / (2 * besselk(p, sqrt(a * b))))
        η = a / 2
        ω = b / 2
        lη = log(η)
        lω = log(ω)
        lx = log(x)
        # calculate first term ρ
        ρ = 0.0
        converged = false
        j = 0
        while !converged && j < 100
            ρ_old = ρ
            ρ += c * (-1)^j * gamma(p - j) * exp((-p + j) * lη + j * lω - lfact(j))
            converged = abs(ρ - ρ_old) < eps()
            j += 1
        end
        # calculate second term σ
        σ = 0.0
        converged = false
        i = 0
        while !converged && i < 100
            σ_old = σ
            j = 0
            k = i
            while j <= i
                l = k * lη + j * lω + (k - j + p) * lx - lfact(k) - lfact(j) 
                σ += (c * (-1)^(k + j + 1) * exp(l)) / (k - j + p)
                j += 1
                k -= 1
            end
            converged = abs(σ - σ_old) < eps()
            i += 1
        end
        1 - (ρ + σ)
    else
        zero(T)
    end
end

function mgf(d::GeneralizedInverseGaussian{T}, t::Real) where {T <: Real}
    if t == zero(t)
        one(T)
    else
        a, b, p = params(d)
        (a / (a - 2t))^(p / 2) * besselk(p, sqrt(b * (a - 2t))) / besselk(p, sqrt(a * b))
    end
end

function cf(d::GeneralizedInverseGaussian{T}, t::Real) where {T <: Real}
    if t == zero(t)
        one(T) + zero(T) * im
    else
        a, b, p = params(d)
        (a / (a - 2t * im))^(p / 2) * besselk(p, sqrt(b * (a - 2t * im))) / besselk(p, sqrt(a * b))
    end
end



#### Sampling

# rand method from:
# Hörmann, W. & J. Leydold. (2014). Generating generalized inverse Gaussian random variates.
# J. Stat. Comput. 24: 547–557. doi:10.1007/s11222-013-9387-3

function rand(d::GeneralizedInverseGaussian)
    a, b, p = params(d)
    α = sqrt(a / b)
    β = sqrt(a * b)
    λ = abs(p)
    if (λ > 1) || (β > 1)
        x = _rou_shift(λ, β)
    elseif β >= min(0.5, (2 / 3) * sqrt(1 - λ))
        x = _rou(λ, β)
    else
        x = _hormann(λ, β)
    end
    p >= 0 ? x / α : 1 / (α * x)
end

function _hormann(λ::Real, β::Real)
    # compute bounding rectangles
    m = β / (1 - λ + sqrt((1 - λ)^2 + β^2))  # mode
    x0 = β / (1 - λ)
    xstar = max(x0, 2 / β)
    # in subdomain (0, x0)
    k1 = exp((λ - 1) * log(m) - β * (m + 1 / m) / 2)
    a1 = k1 * x0
    # in subdomain (x0, 2 / β), may be empty
    if x0 < 2 / β
        k2 = exp(-β)
        a2 = λ == 0 ? k2 * log(2 / (β^2)) : k2 * ((2 / β)^λ - x0^λ) / λ
    else
        k2 = 0
        a2 = 0
    end
    # in subdomain (xstar, Inf)
    k3 = xstar^(λ - 1)
    a3 = 2k3 * exp(-xstar * β / 2) / β
    a = a1 + a2 + a3

    # perform rejection sampling
    while true
        u = rand()
        v = a * rand()
        if v <= a1  # in subdomain (0, x0)
            x = x0 * v / a1
            h = k1
        elseif v <= a1 + a2  # in subdomain (x0, 2 / β)
            v -= a1
            x = λ == 0 ? β * exp(v * exp(β)) : (x0^λ + v * λ / k2)^(1 / λ)
            h = k2 * x^(λ - 1)
        else  # in subdomain (xstar, Inf)
            v -= a1 + a2
            x = -2log(exp(-xstar * β / 2) - v * β / (2k3)) / β
            h = k3 * exp(-x * β / 2)
        end
        if log(u * h) <= ((λ - 1) * log(x) - β * (x + 1 / x) / 2)
            return x
        end
    end
end

function _rou(λ::Real, β::Real)
    # compute bounding rectangle
    m = β / (1 - λ + sqrt((1 - λ)^2 + β^2))  # mode
    k = 0.5(λ - 1) - 0.25β * (m + 1 / m)
    xpos = (1 + λ + sqrt((1 + λ)^2 + β^2)) / β
    upos = xpos * exp(0.5(λ - 1) * log(xpos) - 0.25β * (m + 1 / m) - k)
    
    # perform rejection sampling
    while true
        u = upos * rand()
        v = rand()
        x = u / v
        if log(v) <= (0.25(λ - 1) * log(x) - 0.25β * (x + 1 / x) - k) 
            return x
        end
    end
end

function _rou_shift(λ::Real, β::Real)
    # compute bounding rectangle
    m = (λ - 1 + sqrt((λ - 1)^2 + β^2)) / β  # mode
    a = -2(λ + 1) / β - m
    b = 2(λ - 1) * m / β - 1
    p = b - (a^2) / 3
    q = 2(a^3) / 27 - (a * b) / 3 + m
    ϕ = acos(-(q / (2sqrt(-(p^3) / 27))))  # Cardano's formula
    r = sqrt(-4p / 3)
    xneg = r * cos(ϕ / 3 + 4π / 3) - a / 3
    xpos = r * cos(ϕ / 3) - a / 3
    k = 0.5(λ - 1) * log(m) - 0.25β * (m + 1 / m)
    uneg = (xneg - m) * exp(0.5(λ - 1) * log(xneg) - 0.25β * (xneg + 1 / xneg) - k)
    upos = (xpos - m) * exp(0.5(λ - 1) * log(xpos) - 0.25β * (xpos + 1 / xpos) - k)

    # perform rejection sampling
    while true
        u = (upos - uneg) * rand() + uneg
        v = rand()
        x = u / v + m
        if (x > 0) && (log(v) <= 0.5(λ - 1) * log(x) - 0.25β * (x + 1 / x) - k)
            return x
        end
    end
end

