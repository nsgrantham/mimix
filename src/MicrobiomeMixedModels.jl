module MicrobiomeMixedModels

using Mamba
using Distributions
using YAML
#using RCall

using StatsBase: Histogram, Weights, sample

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
    MIMIXNoFactors,
    fit,
    get_model,
    get_inits,
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

struct MIMIXNoFactors <: MicrobiomeModel end

include("generalized-inverse-gaussian.jl")  # distribution used in MIMIX
include(joinpath("models", "mimix.jl"))
include(joinpath("models", "mimix-no-factors.jl"))
include("utils.jl")

end  # module
