import Turing: ~, sample, MH, NUTS, @model
import Serialization: serialize, deserialize

@model function sfit(measurements = log.(loadMixData(; discard_small = true)."Value"))
    Rtot_dist = importInVitroRtotDist()
    Kav_dist = importKavDist(; inflation = 0.1)
    Kav_dist = Matrix(Kav_dist[:, Not("IgG")])

    Rtot = Vector(undef, length(Rtot_dist))
    Kav = Matrix(undef, size(Kav_dist)...)

    for ii in eachindex(Rtot)
        Rtot[ii] ~ Rtot_dist[ii]
    end
    for ii in eachindex(Kav)
        Kav[ii] ~ Kav_dist[ii]
    end

    f4 ~ f4Dist
    f33 ~ f33Dist
    KxStar ~ KxStarDist

    Rtotd, vals, KxStar, Kav = dismantle_x0(vcat(Rtot, [f4, f33, KxStar], reshape(Kav, :)))
    lsigma = 0.1
    lps = log.(mixturePredictions(; Rtot = Rtotd, Kav = Kav, KxStar = KxStar, vals = vals)."Predict")
    for ii in eachindex(measurements)
        measurements[ii] ~ LogNormal(lps[ii], lsigma)
    end
end

function runMCMC(fname = "MCMC_nuts_1000.dat")
    if isfile(fname)
        return deserialize(fname)
    end
    m = sfit()
    c = sample(m, NUTS(), 1000, discard_initial = 500)
    f = serialize(fname, c)
    return c
end

function plotHistPriorDist(dat, dist, name)
    dat = reshape(dat, :)
    xxs = LinRange(dist.μ - 4 * dist.σ, dist.μ + 4 * dist.σ, 100)
    yys = [pdf(dist, xx) for xx in xxs]
    pl = plot(
        layer(x = xxs, y = yys, Geom.line, color = [colorant"red"], order = 1),
        layer(DataFrame("Value" => dat), x = "Value", Geom.histogram(bincount = 8, density = true)),
        Scale.x_continuous(labels = x -> @sprintf("%.3E", exp(x))),
        Guide.xticks(orientation = :horizontal),
        Guide.xlabel("Value"),
        Guide.ylabel(nothing),
        Guide.title(name),
    )
    return pl
end

function plot_MCMC_dists(c = runMCMC())
    setGadflyTheme()

    # Plot Kav's
    ligg = length(humanIgG)
    lfcr = length(humanFcgRiv)
    Kav_dist = importKavDist(; inflation = 0.1, retdf = false)
    Kav_pls = Matrix{Plot}(undef, ligg, lfcr)
    for ii in eachindex(Kav_pls)
        IgGname = humanIgG[(ii - 1) % ligg + 1]
        FcRname = humanFcgRiv[(ii - 1) ÷ ligg + 1]
        name = IgGname * " to " * FcRname
        Kav_pls[ii] = plotHistPriorDist(c["lKav[$ii]"].data, Kav_dist[ii], name)
    end
    Kav_plot = plotGrid((ligg, lfcr), permutedims(Kav_pls, (2, 1)); sublabels = false)
    draw(SVG("MCMC_Kav.svg", 16inch, 13inch), Kav_plot)

    # Plot Rtot's
    Rtot_pls = Vector{Plot}(undef, lfcr)
    Rtot_dist = importInVitroRtotDist()
    for ii in eachindex(Rtot_pls)
        FcRname = humanFcgRiv[ii]
        Rtot_pls[ii] = plotHistPriorDist(c["lRtot[$ii]"].data, Rtot_dist[ii], FcRname)
    end
    Rtot_plot = plotGrid((1, lfcr), Rtot_pls; sublabels = false)
    draw(SVG("MCMC_Rtot.svg", 16inch, 4inch), Rtot_plot)

    # Plot f4, f33, KxStar
    other_pls = Vector{Plot}(undef, 3)
    other_pls[1] = plotHistPriorDist(c["lf4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plotHistPriorDist(c["lf33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plotHistPriorDist(c["lKxStar"].data, KxStarDist, "K<sub>x</sub><sup>*</sup>")
    other_plot = plotGrid((1, 3), other_pls; sublabels = false)
    draw(SVG("MCMC_others.svg", 8inch, 4inch), other_plot)
end


function MCMCresults2fit(c = runMCMC())
    m = DataFrame(mean(c))
    m."order" .= 0
    m[startswith.(String.(m."parameters"), "lRtot"), "order"] .= 1
    m[startswith.(String.(m."parameters"), "lf"), "order"] .= 2
    m[startswith.(String.(m."parameters"), "lKxStar"), "order"] .= 3
    m[startswith.(String.(m."parameters"), "lKav"), "order"] .= 4
    Rtot, vals, KxStar, Kav = dismantle_x0(exp.(sort(m, "order")."mean"))

    return mixturePredictions(; Rtot = Rtot, Kav = Kav, KxStar = KxStar, vals = vals)
end
