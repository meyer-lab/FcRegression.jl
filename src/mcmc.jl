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

    lf4 ~ Normal(log(4), 0.2*log(4))
    lf33 ~ Normal(log(33), 0.2*log(33))
    lKxStar ~ Normal(log(FcRegression.KxConst), 4.37)
    
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

function plot_MCMC_dists(c::Chains)
    # Kav, Rtot, f4, f33, KxStar
    # StatsPlots
    #histogram(c2["lKav[1]"])

    # Gadfly
    ligg = length(FcRegression.humanIgG)
    lfcr = length(FcRegression.humanFcgRiv)
    Kav_dist = FcRegression.importKavDist(; inflation = 0.1, retdf = false)
    Kav_pls = Matrix{Plot}(undef, ligg, lfcr)
    for ii in eachindex(Kav_pls)
        dat = exp.(reshape(c["lKav[$ii]"].data, :))
        IgGname = FcRegression.humanIgG[(ii-1)%ligg+1]
        FcRname = FcRegression.humanFcgRiv[(ii-1)÷ligg+1]
        name = IgGname * " to " * FcRname
        dist = Kav_dist[ii]
        xxs = (dist.μ - 4 * dist.σ):0.01:(dist.μ + 4 * dist.σ)
        yys = [pdf(dist, xx) * length(dat) for xx in xxs]
        Kav_pls[ii] = plot(
            #layer(x=exp.(xxs), y=yys, Geom.line),
            DataFrame("Value" => dat), x="Value", Geom.histogram(bincount = 8, density = true),
            xintercept=[exp(dist.μ)], Geom.vline,
            Scale.y_continuous(labels = n -> string(n/length(dat))),
            #layer(DataFrame("Value" => dat), x="Value", Geom.histogram(bincount = 8, density=true)),
            Guide.xticks(orientation=:horizontal),
            Guide.title(name)
        )
    end
    Kav_plot = FcRegression.plotGrid((ligg, lfcr), permutedims(Kav_pls, (2 ,1)); sublabels = false);
    draw(SVG("MCMC_Kav.svg", 16inch, 13inch), Kav_plot)
end