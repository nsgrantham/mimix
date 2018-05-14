module MicrobiomeMixedModels

using Mamba
using Distributions
using YAML
#using RCall

using StatsBase: Histogram

import Base: convert
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
    fit

#abstract type MicrobiomeModel end
#struct MIMIX <: MicrobiomeModel end
#struct MIMIXNoFactors <: MicrobiomeModel end

include("generalized-inverse-gaussian.jl")  # distribution used in MIMIX
include(joinpath("models", "mimix.jl"))
include(joinpath("models", "mimix-no-factors.jl"))
include("utils.jl")

end  # module
