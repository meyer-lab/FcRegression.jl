function dist_violin_plot(df, dist_list; title = "", x_name = "Receptor")
    setGadflyTheme()
    
    odf = DataFrame([names(df)[i] => rand(dist_list[i], 500) for i = 1:size(df)[2] if dist_list[i].μ > 0])
    df = stack(df, variable_name=x_name, value_name="Value")
    odf = stack(odf, variable_name=x_name, value_name="Value")

    return plot(layer(df, x=x_name, y="Value", Geom.violin, Theme(default_color=colorant"red")),
        layer(odf, x=x_name, y="Value", Geom.violin, Theme(default_color=colorant"green")),
        Coord.cartesian(ymin=4, ymax=8),
        Scale.y_log10, Guide.manual_color_key("Legend", ["Prior", "Posterior"], ["green", "red"]),
        Guide.ylabel("<i>K</i><sub>a</sub> (M<sup>-1</sup>)"), Guide.title(title))
end

function plot_MCMC_affinity(c = runMCMC())
    df = DataFrame(c)
    df_Rtot = df[500:1000, [contains(n, "Rtot") for n in names(df)]]
    rename!(df_Rtot, humanFcgRiv)
    pl_Rtot = dist_violin_plot(df_Rtot, importInVitroRtotDist(); title = "Distributions of Receptor Expression", x_name="Receptor")

    df_Kav = df[500:1000, [contains(n, "Kav") for n in names(df)]]
    Kav_priors = importKavDist()
    Kav_posts = deepcopy(Kav_priors)
    Kav_posts[!, Not("IgG")] = reshape([df_Kav[!, "Kav[$i]"] for i = 1:size(df_Kav)[2]], size(Kav_priors)[1], :)

    pl_igg = Vector{Union{Gadfly.Plot, Context}}(undef, 0)
    for igg in Kav_priors[!, "IgG"]
        priors = reshape(Matrix(Kav_priors[Kav_priors."IgG" .== igg, Not("IgG")]), :)
        posts = DataFrame(hcat(reshape(Matrix(Kav_posts[Kav_posts."IgG" .== "IgG1", Not("IgG")]), :)...), humanFcgRiv)
        append!(pl_igg, [dist_violin_plot(posts, priors; title = "h$igg Affinities Distributions")])
    end
    
    misc_priors = [f4Dist, f33Dist, f4conv_dist, f33conv_dist]
    df_misc = df[500:1000, ["f4", "f33", "f4conv", "f33conv"]]
    rename!(df_misc, ["Valency\nf = 4", "Valency\nf = 33", "Conv. factor\n(f = 4)", "Conv. factor\n(f = 33)"])
    pl_misc = dist_violin_plot(df_misc, misc_priors; title = "Distributions of Valencies and Conversion factors")

    pl_KxStar = dist_violin_plot(df[500:1000, ["KxStar"]], [KxStarDist]; title = "Distributions of KxStar")
    return pl_Rtot, pl_igg, pl_misc, pl_KxStar
end
