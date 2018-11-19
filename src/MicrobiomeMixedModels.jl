module MicrobiomeMixedModels

using Mamba
using Distributions
using LinearAlgebra
using YAML
using RCall

using StatsBase: Histogram, Weights, sample
using DataFrames: DataFrame, names!

import Base: convert, maximum
import Distributions: params, partype, mean, var, mode, pdf, logpdf, rand
import StatsBase: fit

export
    GeneralizedInverseGaussian,
    params,
    partype,
    mean,
    var,
    mode,
    pdf,
    logpdf,
    rand,
    MIMIX,
    MIMIXGaussian,
    MIMIXNoFactors,
    permanova,
    fit,
    get_model,
    get_inits,
    get_post,
    clr,
    iclr,
    parse_config,
    load_config

abstract type MicrobiomeModel end

struct MIMIX <: MicrobiomeModel
    factors::Int
    function MIMIX(factors::Int)
        @assert factors > zero(factors)
        new(factors)
    end
end

struct MIMIXGaussian <: MicrobiomeModel
    factors::Int
    function MIMIXGaussian(factors::Int)
        @assert factors > zero(factors)
        new(factors)
    end
end

struct MIMIXNoFactors <: MicrobiomeModel end

include("generalized-inverse-gaussian.jl")  # distribution used in MIMIX
include(joinpath("models", "mimix.jl"))
include(joinpath("models", "mimix-gaussian.jl"))
include(joinpath("models", "mimix-no-factors.jl"))
include(joinpath("models", "permanova.jl"))
include("utils.jl")

end  # module
