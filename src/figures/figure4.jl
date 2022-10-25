""" Figure 4: Affinity updates """


""" Making one single violin plot. Called by plotAffinityViolin() """
function plot_distribution_violins(
    df::AbstractDataFrame,
    dist_list::Vector{T};
    title = "",
    x_name = "Receptor",
    y_range = (5, 8),
    legend = true,
) where {T <: Distribution}
    setGadflyTheme()
    @assert size(df)[2] == length(dist_list)

    odf = DataFrame([names(df)[i] => rand(dist_list[i], 100_000) for i = 1:size(df)[2] if dist_list[i].μ > 0])
    df = stack(df, variable_name = x_name, value_name = "Value")
    odf = stack(odf, variable_name = x_name, value_name = "Value")
    df."Receptor" = replace.(df."Receptor", "FcgR" => "FcγR")
    odf."Receptor" = replace.(odf."Receptor", "FcgR" => "FcγR")

    return plot(
        layer(df, x = x_name, y = "Value", Geom.violin, Theme(default_color = colorant"firebrick4")),
        layer(odf, x = x_name, y = "Value", Geom.violin, Theme(default_color = colorant"navajowhite2")),
        Coord.cartesian(ymin = y_range[1], ymax = y_range[2]),
        Scale.y_log10,
        legend ? Guide.manual_color_key("Legend", ["Prior", "Posterior"], ["navajowhite2", "firebrick4"]) : style(key_position = :none),
        Guide.ylabel("<i>K</i><sub>a</sub> (M<sup>-1</sup>)"),
        Guide.title(title),
    )
end

""" Make violin plot of affinities from a chain """
function plotAffinityViolin(c::Chains; murine::Bool, y_range = (4, 9))
    pref = murine ? "m" : "h"
    Kav_priors = importKavDist(; murine = murine, retdf = true)

    len = length(Matrix(Kav_priors[!, Not("IgG")]))
    Kav_posts = deepcopy(Kav_priors)
    Kav = [c["Kav[$i]"].data for i = 1:len]
    Kav_posts[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kav_posts)[1], :))

    pls = Vector{Union{Gadfly.Plot, Context}}(undef, size(Kav_priors)[1])
    for (ii, igg) in enumerate(Kav_priors[!, "IgG"])
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(
            hcat([reshape(Kav_posts[Kav_posts."IgG" .== igg, Not("IgG")][1, i], :) for i = 1:(size(Kav_posts)[2] - 1)]...),
            names(Kav_posts)[2:end],
        )
        pls[ii] = plot_distribution_violins(
            posts,
            priors;
            y_range = y_range,
            title = "$pref$igg Affinity Distributions",
            legend = (ii == length(Kav_priors[!, "IgG"])),
        )
    end
    return pls
end

function figure4(; kwargs...)
    c = rungMCMC("humanKavfit_0701.dat"; dat = :hCHO, mcmc_iter = 1_000)
    pl_igg = plotAffinityViolin(c; murine = false, y_range = (5, 8))

    p1, p2 = bindVSaff(extractNewHumanKav(); affinity_name = "Updated")

    draw(PDF("output/figure4.pdf", 10inch, 6inch), pp)
    )
        kwargs...,
        sublabels = "abcde f ",
        widths = [3 3 3 3.8; 3.8 0.4 5 4.1],
        [pl_igg[1], pl_igg[2], pl_igg[3], pl_igg[4], p1, nothing, p2, nothing];
        (2, 4),
    pp = plotGrid(
end
