import Turing: ~, sample, MH, NUTS, @model
import Serialization: serialize, deserialize
using Distributions
using LinearAlgebra

const f4Dist = LogNormal(log(4), 0.1)
const f33Dist = LogNormal(log(33), 0.1)
const KxStarDist = LogNormal(log(KxConst), 2.0)   # ~ 4.37 in Robinett
const f4conv_dist = LogNormal(log(2.27), 0.2)   # std ~= 0.473
const f33conv_dist = LogNormal(log(3.26), 0.2)  #std ~= 0.672

@model function sfit(df, values; robinett = false)
    Rtot_dist = importInVitroRtotDist(robinett)
    Kav_dist = importKavDist(; inflation = 0.0)
    Kav_dist = Matrix(Kav_dist[:, Not("IgG")])

    Rtot = Vector(undef, length(Rtot_dist))
    Kav = Matrix(undef, size(Kav_dist)...)

    # Order of distribution definitions here matches MAPLikelihood()
    for ii in eachindex(Rtot)
        Rtot[ii] ~ truncated(Rtot_dist[ii], 10, 1E7)
    end
    Rtotd = Dict([humanFcgRiv[ii] => Rtot[ii] for ii = 1:length(humanFcgRiv)])

    for ii in eachindex(Kav)
        Kav[ii] ~ truncated(Kav_dist[ii], 10, 1E9)
    end

    Kavd = deepcopy(importKav(; murine = false, invitro = true, retdf = true))
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(Kav)

    f4 ~ truncated(f4Dist, 1.0, 8.0)
    f33 ~ truncated(f33Dist, 8.0, 50.0)
    KxStar ~ truncated(KxStarDist, 1E-16, 1E-9)
    f4conv ~ truncated(f4conv_dist, 1.0, 4.0)
    f33conv ~ truncated(f33conv_dist, 2.0, 6.0)

    if any(Kav .< 0.0) || (f4 < 0.0) || (f33 < 0.0) || (KxStar < 0.0)
        df = deepcopy(df)
        df."Predict" .= -1000.0
    else
        df = mixturePredictions(deepcopy(df); Rtot = Rtotd, Kav = Kavd, KxStar = KxStar, vals = [f4, f33], convs = [f4conv, f33conv])
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function runMCMC(fname = "MCMC_nuts_1000.dat")
    if isfile(fname)
        return deserialize(fname)
    end
    df = loadMixData()

    # only use single IgG
    df = df[(df."%_1" .== 1.0) .| (df."%_2" .== 1.0), :]
    m = sfit(df, df."Value")
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
    draw(PDF("MCMC_Kav.pdf", 16inch, 13inch), Kav_plot)

    # Plot Rtot's
    Rtot_pls = Vector{Plot}(undef, lfcr)
    Rtot_dist = importInVitroRtotDist()
    for ii in eachindex(Rtot_pls)
        FcRname = humanFcgRiv[ii]
        Rtot_pls[ii] = plotHistPriorDist(c["Rtot[$ii]"].data, Rtot_dist[ii], FcRname)
    end
    Rtot_plot = plotGrid((1, lfcr), Rtot_pls; sublabels = false)
    draw(SVG("MCMC_Rtot.svg", 16inch, 4inch), Rtot_plot)
    draw(PDF("MCMC_Rtot.pdf", 16inch, 4inch), Rtot_plot)

    # Plot f4, f33, KxStar
    other_pls = Vector{Plot}(undef, 3)
    other_pls[1] = plotHistPriorDist(c["f4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plotHistPriorDist(c["f33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plotHistPriorDist(c["KxStar"].data, KxStarDist, "K<sub>x</sub><sup>*</sup>")
    other_plot = plotGrid((1, 3), other_pls; sublabels = false)
    draw(SVG("MCMC_others.svg", 8inch, 4inch), other_plot)
    draw(PDF("MCMC_others.pdf", 8inch, 4inch), other_plot)
end
