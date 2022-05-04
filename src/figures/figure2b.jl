function dist_violin_plot(df::AbstractDataFrame, dist_list::Vector{T}; 
        title = "", x_name = "Receptor", y_range = (4, 8)) where T <: Distribution
    setGadflyTheme()
    @assert size(df)[2] == length(dist_list)

    odf = DataFrame([names(df)[i] => rand(dist_list[i], 100_000) for i = 1:size(df)[2] if dist_list[i].μ > 0])
    df = stack(df, variable_name = x_name, value_name = "Value")
    odf = stack(odf, variable_name = x_name, value_name = "Value")

    return plot(
        layer(df, x = x_name, y = "Value", Geom.violin, Theme(default_color = colorant"red")),
        layer(odf, x = x_name, y = "Value", Geom.violin, Theme(default_color = colorant"green")),
        Coord.cartesian(ymin = y_range[1], ymax = y_range[2]),
        Scale.y_log10,
        Guide.manual_color_key("Legend", ["Prior", "Posterior"], ["green", "red"]),
        Guide.ylabel("<i>K</i><sub>a</sub> (M<sup>-1</sup>)"),
        Guide.title(title),
    )
end

function plot_MCMC_affinity(c = runMCMC())
    df = DataFrame(c)
    df_Rtot = df[!, [contains(n, "Rtot") for n in names(df)]]
    rename!(df_Rtot, humanFcgRiv)
    pl_Rtot = dist_violin_plot(df_Rtot, importInVitroRtotDist(); title = "Distributions of Receptor Expression", x_name = "Receptor")

    df_Kav = df[!, [contains(n, "Kav") for n in names(df)]]
    Kav_priors = importKavDist()
    Kav_posts = deepcopy(Kav_priors)
    Kav_posts[!, Not("IgG")] = reshape([df_Kav[!, "Kav[$i]"] for i = 1:size(df_Kav)[2]], size(Kav_priors)[1], :)

    pl_igg = Vector{Union{Gadfly.Plot, Context}}(undef, 0)
    for igg in Kav_priors[!, "IgG"]
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(hcat(reshape(Matrix(Kav_posts[Kav_posts."IgG" .== "IgG1", Not("IgG")]), :)...), humanFcgRiv)
        append!(pl_igg, [dist_violin_plot(posts, priors; title = "h$igg Affinities Distributions")])
    end

    misc_priors = [f4Dist, f33Dist, KxStarDist]
    df_misc = df[!, ["f4", "f33", "KxStar"]]
    rename!(df_misc, ["Valency\nf = 4", "Valency\nf = 33", "KxStar"])
    pl_misc = dist_violin_plot(df_misc, misc_priors; title = "Distributions of Valencies and KxStar")
    return pl_Rtot, pl_igg, pl_misc
end
