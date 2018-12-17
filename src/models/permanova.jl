
function permanova(data)
    X = data[:X]
    Y = data[:Y]
    Z = data[:Z]
    R"""
    suppressMessages(suppressWarnings(library(vegan)))
    z <- data.frame($Z)
    z[] <- lapply(z, as.factor)
    x <- data.frame($X)
    x[] <- lapply(x, as.factor)
    sim <- adonis($Y ~ $X + $Z, method="bray", nperm=9999)
    print(sim)
    pval <- sim$aov.tab$`Pr(>F)`[1]
    """
    @rget pval
    return Dict{Symbol, Any}(:pval => pval)
end