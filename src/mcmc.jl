import Turing: ~, sample, MH, NUTS, @model
import Serialization: serialize, deserialize
using Distributions
using ForwardDiff
using LinearAlgebra
import Statistics: cor


@model function sfit(df)
    df = deepcopy(df)

    Rtot_dist = importInVitroRtotDist()
    Kav_dist = importKavDist(; inflation = 0.1)
    Kav_dist = Matrix(Kav_dist[:, Not("IgG")])

    Rtot = Vector(undef, length(Rtot_dist))
    Kav = Matrix(undef, size(Kav_dist)...)

    for ii in eachindex(Rtot)
        Rtot[ii] ~ truncated(Rtot_dist[ii], 1000, 1E7)
    end
    for ii in eachindex(Kav)
        Kav[ii] ~ truncated(Kav_dist[ii], 100, 1E9)
    end

    f4 ~ LogNormal(log(4), 0.1)
    f33 ~ LogNormal(log(33), 0.1)
    KxStar ~ truncated(LogNormal(log(KxConst), 0.1), 1E-14, 1E-10)

    x0 = vcat(Rtot, [f4, f33, KxConst], reshape(Kav, :)) # Simplifying fitting for a bit
    T = typeof(x0[1])
    Rtotd, vals, KxStar, Kav = dismantle_x0(T.(x0))
    df = mixturePredictions(df; Rtot = Rtotd, Kav = Kav, KxStar = KxStar, vals = vals)

    lps = df."Predict"
    measurements = df."Value"
    @assert all(isfinite(lps))
    @assert all(isfinite(measurements))
    @assert size(lps) == size(measurements)

    # Least squares with one var and no intercept
    scale = sum(measurements .* lps) / sum(measurements .* measurements)
    diff = log.(lps .* scale) - log.(measurements)

    # println(ForwardDiff.value(cor(log.(lps .* scale), log.(measurements))))
    # Shows a correlation of 0.51
    
    # diffSum = norm(diff / 0.1) # Scale by expected variance
    diff ~ MvNormal(zeros(size(diff)), 1.0 * I)
    # diffSum ~ Chi(length(diff)) # Speeds up a lot
    nothing
end

function runMCMC(fname = "MCMC_nuts_1000.dat")
    if isfile(fname)
        return deserialize(fname)
    end
    df = loadMixData()
    m = sfit(df)
    c = sample(m, NUTS(), 1_000)
    f = serialize(fname, c)
    return c
end

function plotHistPriorDist(dat, dist, name)
    dat = reshape(dat, :)
    xxs = LinRange(dist.μ - 4 * dist.σ, dist.μ + 4 * dist.σ, 100)
    yys = [pdf(dist, exp.(xx)) for xx in xxs]
    yys = yys ./ maximum(yys)
    pl = plot(
        layer(x = xxs, y = yys, Geom.line, color = [colorant"red"], order = 1),
        layer(DataFrame("Value" => log.(dat)), x = "Value", Geom.histogram(bincount = 20, density = true)),
        Scale.x_continuous(labels = x -> @sprintf("%.1E", exp(x))),
        Guide.xticks(orientation = :horizontal),
        Guide.xlabel("Value"),
        Guide.ylabel(nothing),
        Guide.title(name),
    )
    return pl
end

function plot_MCMC_dists(c = runMCMC())
    setGadflyTheme()
    c = c[500:1000]

    # Plot Kav's
    ligg = length(humanIgG)
    lfcr = length(humanFcgRiv)
    Kav_dist = importKavDist(; inflation = 0.1, retdf = false)
    Kav_pls = Matrix{Plot}(undef, ligg, lfcr)
    for ii in eachindex(Kav_pls)
        IgGname = humanIgG[(ii - 1) % ligg + 1]
        FcRname = humanFcgRiv[(ii - 1) ÷ ligg + 1]
        name = IgGname * " to " * FcRname
        Kav_pls[ii] = plotHistPriorDist(c["Kav[$ii]"].data, Kav_dist[ii], name)
    end
    Kav_plot = plotGrid((ligg, lfcr), permutedims(Kav_pls, (2, 1)); sublabels = false)
    draw(SVG("MCMC_Kav.svg", 16inch, 13inch), Kav_plot)

    # Plot Rtot's
    Rtot_pls = Vector{Plot}(undef, lfcr)
    Rtot_dist = importInVitroRtotDist()
    for ii in eachindex(Rtot_pls)
        FcRname = humanFcgRiv[ii]
        Rtot_pls[ii] = plotHistPriorDist(c["Rtot[$ii]"].data, Rtot_dist[ii], FcRname)
    end
    Rtot_plot = plotGrid((1, lfcr), Rtot_pls; sublabels = false)
    draw(SVG("MCMC_Rtot.svg", 16inch, 4inch), Rtot_plot)

    # Plot f4, f33, KxStar
    other_pls = Vector{Plot}(undef, 3)
    other_pls[1] = plotHistPriorDist(c["f4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plotHistPriorDist(c["f33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plotHistPriorDist(c["KxStar"].data, KxStarDist, "K<sub>x</sub><sup>*</sup>")
    other_plot = plotGrid((1, 3), other_pls; sublabels = false)
    draw(SVG("MCMC_others.svg", 8inch, 4inch), other_plot)
end
