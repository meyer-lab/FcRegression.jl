import Turing: ~, sample, MH, NUTS, @model
import Serialization: serialize, deserialize
using Distributions
using LinearAlgebra

const f4Dist = LogNormal(log(4), 0.2)   # std ~= 0.82
const f33Dist = LogNormal(log(33), 0.2)   # std ~= 6.80
const KxStarDist = LogNormal(log(KxConst), 2.0)   # std ~= 4.37 in Robinett
const f4conv_dist = LogNormal(log(2.27), 0.2)   # std ~= 0.468
const f33conv_dist = LogNormal(log(3.26), 0.2)  #std ~= 0.672

@model function sfit(df, values; robinett = false, Kavd = importKav(; murine = false, invitro = true, retdf = true))
    Rtot_dist = importInVitroRtotDist(robinett)
    Kav_dist = importKavDist()
    Kav_dist = Matrix(Kav_dist[:, Not("IgG")])

    Rtot = Vector(undef, length(Rtot_dist))
    Kav = Matrix(undef, size(Kav_dist)...)

    # Order of distribution definitions here matches MAPLikelihood()
    for ii in eachindex(Rtot)
        Rtot[ii] ~ truncated(Rtot_dist[ii], 10, 1E8)
    end
    Rtotd = Dict([humanFcgRiv[ii] => Rtot[ii] for ii = 1:length(humanFcgRiv)])

    for ii in eachindex(Kav)
        Kav[ii] ~ truncated(Kav_dist[ii], 10, 1E10)
    end

    if !robinett    # Don't fit affinity for Robinett data
        Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(Kav)
    end

    f4 ~ truncated(f4Dist, 1.0, 8.0)
    f33 ~ truncated(f33Dist, 8.0, 50.0)
    KxStar ~ truncated(KxStarDist, 1E-18, 1E-9)
    f4conv ~ truncated(f4conv_dist, 1.0, 5.0)
    f33conv ~ truncated(f33conv_dist, 2.0, 10.0)

    if all(Rtot .> 0.0) && all(Kav .> 0.0) && all([f4, f33, KxStar, f4conv, f33conv] .> 0.0)
        df = mixturePredictions(deepcopy(df); Rtot = Rtotd, Kav = Kavd, KxStar = KxStar, vals = [f4, f33], convs = [f4conv, f33conv])
    else
        df = deepcopy(df)
        df."Predict" .= -1000.0
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function runMCMC(fname = "MCMC_nuts_wconvs_0405.dat")
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

""" Making a single subplot for priors and posteriors """
function plotHistPriorDist(dat, dist, name)
    dat = reshape(dat, :)
    xxs = exp.(LinRange(dist.μ - 4 * dist.σ, dist.μ + 4 * dist.σ, 100))
    yys = [pdf(dist, xx) for xx in xxs]
    yys = yys ./ maximum(yys)
    pl = plot(
        layer(x = xxs, y = yys, Geom.line, color = [colorant"red"], order = 1),
        layer(DataFrame("Value" => dat), x = "Value", Geom.histogram(bincount = 20, density = true)),
        Scale.x_log10(),
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
    Kav_dist = importKavDist(; retdf = false)
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
    other_pls = Vector{Plot}(undef, 5)
    other_pls[1] = plotHistPriorDist(c["f4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plotHistPriorDist(c["f33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plotHistPriorDist(c["KxStar"].data, KxStarDist, "K<sub>x</sub><sup>*</sup>")
    other_pls[4] = plotHistPriorDist(c["f4conv"].data, f4conv_dist, "f = 4 conversion factor")
    other_pls[5] = plotHistPriorDist(c["f33conv"].data, f33conv_dist, "f = 33 conversion factor")
    other_plot = plotGrid((1, 5), other_pls; sublabels = false)
    draw(SVG("MCMC_others.svg", 12inch, 4inch), other_plot)
    draw(PDF("MCMC_others.pdf", 12inch, 4inch), other_plot)
end

function extractMCMCresults(c = runMCMC("MCMC_nuts_wconvs_0328.dat"))
    c = c[500:1000]
    Rtot = [median(c["Rtot[$i]"].data) for i = 1:length(humanFcgRiv)]
    Rtotd = Dict([humanFcgRiv[ii] => Rtot[ii] for ii = 1:length(humanFcgRiv)])

    Kavd = importKav(; murine = false, invitro = true, retdf = true)
    Kav = [median(c["Kav[$i]"].data) for i = 1:length(importKav(; murine = false, invitro = true, retdf = false))]
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
    
    f4 = median(c["f4"].data)
    f33 = median(c["f33"].data)
    KxStar = median(c["KxStar"].data)
    f4conv = median(c["f4conv"].data)
    f33conv = median(c["f33conv"].data)

    return Rtotd, Kavd, [f4, f33, KxStar, f4conv, f33conv]
end

function MCMC_params_predict(c = runMCMC(), df = loadMixData())
    Rtotd, Kavd, (f4, f33, KxStar, f4conv, f33conv) = extractMCMCresults(c)

    if !("xmin" in names(df))
        df = averageMixData(df)
    end
    return mixturePredictions(df; Rtot = Rtotd, Kav = Kavd, KxStar = KxStar, vals = [f4, f33], convs = [f4conv, f33conv])
end
