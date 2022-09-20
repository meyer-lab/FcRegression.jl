""" A general plotting function for Adjusted vs. Predicted plots """
function plotPredvsMeasured(
    df;
    xx = "Value",
    yy = "Predict",
    xxlabel = "Actual",
    yylabel = "Predicted",
    color = "Receptor",
    shape = "Valency",
    title = "Predicted vs Actual",
    R2pos = (3, 1),
)
    setGadflyTheme()
    df = deepcopy(df)
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    # Move 0 to a small nonzero only when plotting
    df[(df[!, xx]) .<= 0.0, xx] .= minimum(df[(df[!, xx]) .> 0.0, xx]) / 10
    df[(df[!, yy]) .<= 0.0, yy] .= minimum(df[(df[!, yy]) .> 0.0, yy]) / 10
    r2 = R2((df[!, xx]), (df[!, yy]))
    println(r2)
    errbar = "xmin" in names(df)
    return plot(
        df,
        x = xx,
        y = yy,
        xmin = (errbar ? "xmin" : xx),
        xmax = (errbar ? "xmax" : xx),
        color = color,
        shape = shape,
        Geom.point,
        "xmin" in names(df) ? Geom.errorbar : Guide.xlabel(xxlabel),
        Guide.ylabel(yylabel, orientation = :vertical),
        Guide.title(title),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(
            Scale.color_discrete().f(10)[1],
            Scale.color_discrete().f(10)[3],
            Scale.color_discrete().f(10)[2],
            Scale.color_discrete().f(10)[4:end]...,
        ),
        Geom.abline(color = "black"),
        Guide.annotation(
            compose(
                context(),
                text(R2pos[1], R2pos[2], "<i>R</i><sup>2</sup> = " * @sprintf("%.4f", r2)),
                stroke("black"),
                fill("black"),
                font("Helvetica-Bold"),
            ),
        ),
        style(errorbar_cap_length = 0px),
    )
end


""" Making a single subplot for priors and posteriors. Called by plotMCMCdists() """
function plot_dist_histogram(dat::Array{Float64}, dist::Distribution, name::String = ""; bincount = 20)
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


""" Plot prior and posterior distributions from a chain. Will generate three figure files """
function plotMCMCdists(c::Chains, fname::String = ""; murine::Bool)
    setGadflyTheme()
    pnames = [String(s) for s in c.name_map[1]]

    local Kav_dist, IgGs
    Kav_dist = importKavDist(; murine = murine, retdf = false)
    IgGs = importKavDist(; murine = murine, retdf = true)."IgG"
    FcgRs = murine ? murineFcgR : humanFcgRiv
    pref = murine ? "m" : "h"
    ligg, lfcr = length(IgGs), length(FcgRs)

    # Plot Kav's
    if "Kav[1]" in pnames
        Kav_pls = Matrix{Plot}(undef, ligg, lfcr)
        for ii in eachindex(Kav_pls)
            IgGname = IgGs[(ii - 1) % ligg + 1]
            FcRname = FcgRs[(ii - 1) ÷ ligg + 1]
            name = "$pref$IgGname to $pref$FcRname"
            Kav_pls[ii] = plot_dist_histogram(c["Kav[$ii]"].data, Kav_dist[ii], name)
        end
        Kav_plot = plotGrid((ligg, lfcr), permutedims(Kav_pls, (2, 1)); sublabels = false)
        draw(PDF("MCMC_Kav_$fname.pdf", 11inch, 8inch), Kav_plot)
    end

    # Plot Rtot's
    ## TODO: simplify this
    if murine
        cellTypes = unique(importMurineLeukocyte(; average = true)."ImCell")
        lcell = length(cellTypes)
        Rtotd = importRtotDist(:mLeuk; retdf = true)
        Rtotd = Matrix(Rtotd[!, names(Rtotd)[in(cellTypes).(names(Rtotd))]])
        Rtot_pls = Matrix{Plot}(undef, lcell, lfcr)
        for ii in eachindex(Rtot_pls)
            cellname = cellTypes[(ii - 1) % lcell + 1]
            FcRname = murineFcgR[(ii - 1) ÷ lcell + 1]
            name = "m$FcRname on $cellname"
            Rtot_pls[ii] = plot_dist_histogram(c["Rtot[$ii]"].data, Rtotd[ii], name)
        end
        Rtot_plot = plotGrid((lcell, lfcr), permutedims(Rtot_pls, (2, 1)); sublabels = false)
        draw(PDF("MCMC_Rtot_$fname.pdf", 11inch, 14inch), Rtot_plot)
    else # human
        Rtot_pls = Vector{Plot}(undef, lfcr)
        Rtot_dist = importRtotDist(:hCHO)
        for ii in eachindex(Rtot_pls)
            FcRname = humanFcgRiv[ii]
            Rtot_pls[ii] = plot_dist_histogram(c["Rtot[$ii]"].data, Rtot_dist[ii], FcRname)
        end
        Rtot_plot = plotGrid((1, lfcr), Rtot_pls; sublabels = false)
        draw(PDF("MCMC_Rtot_$fname.pdf", 16inch, 4inch), Rtot_plot)
    end

    # Plot f4, f33, KxStar
    ## TODO: add case where not all three parameters appear
    other_pls = Vector{Plot}(undef, 3)
    other_pls[1] = plot_dist_histogram(c["f4"].data, f4Dist, "f = 4 effective valency")
    other_pls[2] = plot_dist_histogram(c["f33"].data, f33Dist, "f = 33 effective valency")
    other_pls[3] = plot_dist_histogram(c["KxStar"].data, KxStarDist, "KxStar")
    other_plot = plotGrid((1, 3), other_pls; sublabels = false)
    draw(PDF("MCMC_others_$fname.pdf", 9inch, 3inch), other_plot)
end


function plotMCMCPredict(c, df::AbstractDataFrame; dat::Symbol, Kav::Union{Nothing, AbstractDataFrame} = nothing, kwargs...)
    p = extractMCMC(c; dat = dat)
    Kavd = (Kav !== nothing) ? Kav : p["Kav"]
    if (p["Kav"] !== nothing) && (Kav isa AbstractDataFrame)
        @warn "Using provided Kav, even though Kav's are fitted in MCMC chain"
    end

    if !("xmin" in names(df))
        df = averageMixData(df)
    end
    ndf = predMix(df; Kav = Kavd, KxStar = p["KxStar"], Rtot = p["Rtot"], fs = [p["f4"], p["f33"]])
    return plotPredvsMeasured(
        ndf;
        xx = "Value",
        yy = "Predict",
        color = (("ImCell" in names(df)) ? "ImCell" : "Receptor"),
        shape = (("Subclass" in names(df)) ? "Subclass" : "Valency"),
        kwargs...,
    )
end

function plotMAPPredict(df::AbstractDataFrame; dat::Symbol)
    m = gmodel(df, df."Value"; dat = dat)
    opts = Optim.Options(iterations = 1000, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    return plotMCMCPredict(opt, df; dat = dat, title = "$dat MAP fitting results")
end
