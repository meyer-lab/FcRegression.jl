import Turing: ~, sample, NUTS, @model, Chains, ADVI, vi, rand, MAP, MLE
import Serialization: serialize, deserialize
using Distributions
using LinearAlgebra


const f4Dist = LogNormal(log(4), 0.2)   # std ~= 0.82
const f33Dist = LogNormal(log(33), 0.2)   # std ~= 6.80
const KxStarDist = LogNormal(log(KxConst), 2.0)   # std ~= 4.37 in Robinett

@model function sfit(df, values; robinett = false, Kavd = importKav(; murine = false, invitro = true, retdf = true))
    Rtot_dist = importInVitroRtotDist(robinett)
    Kav_dist = importKavDist()
    Kav_dist = Matrix(Kav_dist[:, Not("IgG")])

    Rtot = Vector(undef, length(Rtot_dist))
    Kav = Matrix(undef, size(Kav_dist)...)

    # Order of distribution definitions here matches MAPLikelihood()
    for ii in eachindex(Rtot)
        Rtot[ii] ~ Rtot_dist[ii]
    end
    Rtotd = Dict([humanFcgRiv[ii] => Rtot[ii] for ii = 1:length(humanFcgRiv)])

    for ii in eachindex(Kav)
        Kav[ii] ~ Kav_dist[ii]
    end

    if !robinett    # Don't fit affinity for Robinett data
        Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(Kav)
    end

    f4 ~ f4Dist
    f33 ~ f33Dist
    KxStar ~ KxStarDist

    # fit predictions
    if all(Rtot .>= 0.0) && all(0.0 .<= Kav .< Inf) && all(0.0 .< [f4, f33, KxStar] .< Inf)
        df = predictMix(deepcopy(df); recepExp = Rtotd, Kav = Kavd, KxStar = KxStar, vals = [f4, f33])
    else
        df = deepcopy(df)
        df."Predict" .= Inf
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function runMCMC(fname = "humanNUTSfit_0505.dat"; mcmc_iter = 1_000)
    if isfile(fname)
        return deserialize(fname)
    end
    df = loadMixData()

    # only use single IgG
    df = df[(df."%_1" .== 1.0) .| (df."%_2" .== 1.0), :]
    m = sfit(df, df."Value")

    opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)

    f = serialize(fname, c)
    return c
end

function MAPLikelihood(df; robinett = false)
    model = sfit(df, df."Value"; robinett = robinett)
    opts = Optim.Options(iterations = 1000, show_every = 10, show_trace = true)

    opt = optimize(model, MAP(), LBFGS(; m = 20), opts)
    x = opt.values.array

    Rtot = Dict([humanFcgRiv[ii] => x[ii] for ii = 1:length(humanFcgRiv)])
    Kav = deepcopy(importKav(; murine = false, invitro = true, retdf = true))
    Kav[!, Not("IgG")] = reshape(x[10:33], length(humanIgG), length(humanFcgRiv))

    return predictMix(df; recepExp = Rtot, Kav = Kav, KxStar = x[9], vals = [x[7], x[8]])
end

""" Making a single subplot for priors and posteriors """
function plotHistPriorDist(dat::Array{Float64}, dist::Distribution, name::String = ""; bincount = 20)
    dat = reshape(dat, :)
    xxs = exp.(LinRange(dist.μ - 4 * dist.σ, dist.μ + 4 * dist.σ, 100))
    yys = [pdf(dist, xx) for xx in xxs]
    yys = yys ./ maximum(yys)
    pl = plot(
        layer(x = xxs, y = yys, Geom.line, color = [colorant"red"], order = 1),
        layer(DataFrame("Value" => dat), x = "Value", Geom.histogram(bincount = bincount, density = true)),
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
    draw(PDF("MCMC_Kav.pdf", 16inch, 13inch), Kav_plot)

    # Plot Rtot's
    Rtot_pls = Vector{Plot}(undef, lfcr)
    Rtot_dist = importInVitroRtotDist()
    for ii in eachindex(Rtot_pls)
        FcRname = humanFcgRiv[ii]
        Rtot_pls[ii] = plotHistPriorDist(c["Rtot[$ii]"].data, Rtot_dist[ii], FcRname)
    end
    Rtot_plot = plotGrid((1, lfcr), Rtot_pls; sublabels = false)
    draw(PDF("MCMC_Rtot.pdf", 16inch, 4inch), Rtot_plot)

    # Plot f4, f33, KxStar
    other_pls = Vector{Plot}(undef, 3)
    other_pls[1] = plotHistPriorDist(c["f4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plotHistPriorDist(c["f33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plotHistPriorDist(c["KxStar"].data, KxStarDist, "K<sub>x</sub><sup>*</sup>")
    other_plot = plotGrid((1, 3), other_pls; sublabels = false)
    draw(PDF("MCMC_others.pdf", 12inch, 4inch), other_plot)
end

function extractMCMCresults(c = runMCMC())
    Rtot = [median(c["Rtot[$i]"].data) for i = 1:length(humanFcgRiv)]
    Rtotd = Dict([humanFcgRiv[ii] => Rtot[ii] for ii = 1:length(humanFcgRiv)])

    Kavd = importKav(; murine = false, invitro = true, retdf = true)
    Kav = [median(c["Kav[$i]"].data) for i = 1:length(importKav(; murine = false, invitro = true, retdf = false))]
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))

    f4 = median(c["f4"].data)
    f33 = median(c["f33"].data)
    KxStar = median(c["KxStar"].data)

    return Rtotd, Kavd, [f4, f33, KxStar]
end

function MCMC_params_predict_plot(c = runMCMC(), df = loadMixData(); kwargs...)
    Rtotd, Kavd, (f4, f33, KxStar) = extractMCMCresults(c)

    if !("xmin" in names(df))
        df = averageMixData(df)
    end
    ndf = predictMix(df; recepExp = Rtotd, Kav = Kavd, KxStar = KxStar, vals = [f4, f33])
    return plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", color = "Cell", shape = "Valency", title = "Human NUTS fitting results", kwargs...)
end
