""" Making one single violin plot. Called by plotAffinityViolin() """
function plot_distribution_violins(df::AbstractDataFrame, dist_list::Vector{T}; title = "", x_name = "Receptor", y_range = (5, 8)) where {T <: Distribution}
    setGadflyTheme()
    @assert size(df)[2] == length(dist_list)

    odf = DataFrame([names(df)[i] => rand(dist_list[i], 100_000) for i = 1:size(df)[2] if dist_list[i].μ > 0])
    df = stack(df, variable_name = x_name, value_name = "Value")
    odf = stack(odf, variable_name = x_name, value_name = "Value")

    return plot(
        layer(df, x = x_name, y = "Value", Geom.violin, Theme(default_color = colorant"firebrick4")),
        layer(odf, x = x_name, y = "Value", Geom.violin, Theme(default_color = colorant"navajowhite2")),
        Coord.cartesian(ymin = y_range[1], ymax = y_range[2]),
        Scale.y_log10,
        Guide.manual_color_key("Legend", ["Prior", "Posterior"], ["navajowhite2", "firebrick4"]),
        Guide.ylabel("<i>K</i><sub>a</sub> (M<sup>-1</sup>)"),
        Guide.title(title),
    )
end

""" Make violin plot of affinities from a chain """
function plotAffinityViolin(c::Chains; murine::Bool)
    pref = murine ? "m" : "h"
    y_range = murine ? (4, 9) : (5, 8)
    Kav_priors = importKavDist(; murine = murine, retdf = true)

    len = length(Matrix(Kav_priors[!, Not("IgG")]))
    Kav_posts = deepcopy(Kav_priors)
    Kav = [c["Kav[$i]"].data for i = 1:len]
    Kav_posts[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kav_posts)[1], :))
    
    pls = Vector{Union{Gadfly.Plot, Context}}(undef, size(Kav_priors)[1])
    for (ii, igg) in enumerate(Kav_priors[!, "IgG"])
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(hcat([reshape(Kav_posts[Kav_posts."IgG" .== igg, Not("IgG")][1, i], :) for i = 1:(size(Kav_posts)[2]-1)]...), names(Kav_posts)[2:end])
        pls[ii] = plot_distribution_violins(posts, priors; y_range = y_range, title = "$pref$igg Affinities Distributions")
    end
    return pls
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
    return plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", color = (("ImCell" in names(df)) ? "ImCell" : "Receptor"), 
        shape = (("Subclass" in names(df)) ? "Subclass" : "Valency"), kwargs...)
end
