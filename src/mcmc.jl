import Turing: ~, sample, NUTS, @model, Chains, rand, MAP, MLE

""" A generic function to extract fitted parameters from an MCMC """
function extractMCMC(c::Chains; murine::Bool)
    pnames = [String(s) for s in c.name_map[1]]
    out = Dict{String, Union{Number, Dict, DataFrame, Nothing}}()
    if "Rtot[1]" in pnames
        local Rtotd
        if murine
            cells = unique(importMurineLeukocyte()."Cell")
            Rtotd = importCellRtotDist()
            Rtotd = Rtotd[!, ["Receptor"; names(Rtotd)[in(cells).(names(Rtotd))]]]
            Rtot = [median(c["Rtot[$i]"].data) for i = 1:24]
            Rtotd[!, Not("Receptor")] = typeof(Rtot[1, 1]).(reshape(Rtot, size(Rtotd)[1], :))
        else # human
            Rtot = [median(c["Rtot[$i]"].data) for i = 1:length(humanFcgRiv)]
            Rtotd = Dict([humanFcgRiv[ii] => Rtot[ii] for ii = 1:length(humanFcgRiv)])
        end
        out["Rtot"] = Rtotd
    end
    if "Kav[1]" in pnames
        if murine
            Kavd = murineKavDist(; regularKav = true, retdf = true)
        else # human
            Kavd = importKav(; murine = false, invitro = true, retdf = true)
        end
        Kav = [median(c["Kav[$i]"].data) for i = 1:sum(startswith.(pnames, "Kav"))]
        Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
        out["Kav"] = Kavd
    else
        out["Kav"] = nothing
    end
    for var in ["f4", "f33", "KxStar"]
        if var in pnames
            out[var] = median(c[var].data)
        end
    end
    return out
end


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
    Kav_priors = murine ? murineKavDist(; retdf = true) : importKavDist(; retdf = true)

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


