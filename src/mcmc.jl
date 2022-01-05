using Turing
using Serialization


@model function sfit(measurements = FcRegression.loadMixData(; discard_small = true)."Value")
    Rtot_dist = FcRegression.importInVitroRtotDist()
    Kav_dist = FcRegression.importKavDist(; inflation = 0.1)
    Kav_dist = Matrix(Kav_dist[:, Not("IgG")])

    lRtot = Vector(undef, length(Rtot_dist))
    lKav = Matrix(undef, size(Kav_dist)...)
    
    
    for ii in eachindex(lRtot)
        lRtot[ii] ~ Rtot_dist[ii]
    end
    for ii in eachindex(lKav)
        lKav[ii] ~ Kav_dist[ii]
    end

    lf4 ~ FcRegression.f4Dist
    lf33 ~ FcRegression.f33Dist
    lKxStar ~ FcRegression.KxStarDist
    
    Rtotd, vals, KxStar, Kav = FcRegression.dismantle_x0(exp.(vcat(lRtot, [lf4, lf33, lKxStar], reshape(lKav, :))))
    sigma = 0.1
    ps = FcRegression.mixturePredictions(;
        Rtot = Rtotd,
        Kav = Kav,
        KxStar = KxStar,
        vals = vals,
    )."Predict"
    for ii in eachindex(measurements)
        measurements[ii] ~ Normal(ps[ii], ps[ii] * sigma)
    end
end

function runMCMC(fname = "MCMC_run_100000.dat")
    if isfile(fname)
        return deserialize(fname)
    end
    m = sfit()
    c = sample(m, MH(), 1000, discard_initial = 500)
    f = serialize(fname, c)
    return c
end

function plotHistPriorDist(dat, dist, name)
    dat = reshape(dat, :)
    xxs = LinRange(dist.μ - 4 * dist.σ, dist.μ + 4 * dist.σ, 100)
    yys = [pdf(dist, xx) for xx in xxs]
    pl = plot(
        layer(x=xxs, y=yys, Geom.line, color=[colorant"red"], order=1),
        layer(DataFrame("Value" => dat), x="Value", Geom.histogram(bincount = 8, density=true)),
        Scale.x_continuous(labels = x -> @sprintf("%.3E", exp(x))),
        Guide.xticks(orientation=:horizontal),
        Guide.
        Guide.title(name)
    )
    return pl
end

function plot_MCMC_dists()
    c = runMCMC()

    # Plot Kav's
    ligg = length(FcRegression.humanIgG)
    lfcr = length(FcRegression.humanFcgRiv)
    Kav_dist = FcRegression.importKavDist(; inflation = 0.1, retdf = false)
    Kav_pls = Matrix{Plot}(undef, ligg, lfcr)
    for ii in eachindex(Kav_pls)
        IgGname = FcRegression.humanIgG[(ii-1)%ligg+1]
        FcRname = FcRegression.humanFcgRiv[(ii-1)÷ligg+1]
        name = IgGname * " to " * FcRname
        Kav_pls[ii] = plotHistPriorDist(c["lKav[$ii]"].data, Kav_dist[ii], name)
    end
    Kav_plot = FcRegression.plotGrid((ligg, lfcr), permutedims(Kav_pls, (2 ,1)); sublabels = false);
    draw(SVG("MCMC_Kav.svg", 16inch, 13inch), Kav_plot)

    # Plot Rtot's
    Rtot_pls = Vector{Plot}(undef, lfcr)
    Rtot_dist = FcRegression.importInVitroRtotDist()
    for ii in eachindex(Rtot_pls)
        FcRname = FcRegression.humanFcgRiv[ii]
        Rtot_pls[ii] = plotHistPriorDist(c["lRtot[$ii]"].data, Rtot_dist[ii], FcRname)
    end
    Rtot_plot = FcRegression.plotGrid((1, lfcr), Rtot_pls; sublabels = false);
    draw(SVG("MCMC_Rtot.svg", 16inch, 4inch), Rtot_plot)

    # Plot f4, f33, KxStar
    other_pls = Vector{Plot}(undef, 3)
    other_pls[1] = plotHistPriorDist(c["lf4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plotHistPriorDist(c["lf33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plotHistPriorDist(c["lKxStar"].data, KxStarDist, "K<sub>x</sub><sup>*</sup>")
    other_plot = FcRegression.plotGrid((1, 3), other_pls; sublabels = false);
    draw(SVG("MCMC_others.svg", 8inch, 4inch), other_plot)
end

