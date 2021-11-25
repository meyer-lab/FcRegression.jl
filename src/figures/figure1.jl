""" Figure 1: show the mixture IC binding data """

using Printf
using ColorSchemes

""" Original measurements with middle 50% as error bar """
function splot_origData(df; match_y = true)
    cell = unique(df."Cell")[1]
    IgGX = unique(df."subclass_1")[1]
    IgGY = unique(df."subclass_2")[1]
    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]

    ymax = Dict(
        "FcgRI" => 8e3,
        "FcgRIIA-131H" => 2.5e4,
        "FcgRIIA-131R" => 2.5e4,
        "FcgRIIB-232I" => 2.5e3,
        "FcgRIIIA-158F" => 2.5e4,
        "FcgRIIIA-158V" => 1e4,
    )
    return plot(
        df,
        x = "%_1",
        y = "Median",
        ymin = "xmin",
        ymax = "xmax",
        color = "Valency",
        Geom.point,
        Geom.line,
        Geom.errorbar,
        Scale.x_continuous(labels = n -> "$IgGX $(n*100)%\n$IgGY $(100-n*100)%"),
        Scale.y_continuous(; maxvalue = match_y ? ymax[cell] : maximum(df."xmax")),
        Scale.color_discrete_manual(palette[1], palette[2]),
        Guide.xlabel("", orientation = :horizontal),
        Guide.ylabel("RFU", orientation = :vertical),
        Guide.xticks(orientation = :horizontal),
        Guide.title("$IgGX-$IgGY in $cell"),
    )
end

function bindVSaff()
    hKav = importKav(; murine = false, retdf = true)

    # Binding data, keep single IgG subclass only
    df = averageMixData(loadMixData())
    df = df[(df."%_1" .== 1.0) .| (df."%_2" .== 1.0), :]
    df."Subclass" = [r."%_1" >= 1 ? r."subclass_1" : r."subclass_2" for r in eachrow(df)]
    df = df[!, ["Valency", "Cell", "Subclass", "Value"]]
    lower(x) = quantile(x, 0.25)
    upper(x) = quantile(x, 0.75)
    df = combine(
        groupby(df, ["Valency", "Cell", "Subclass"]),
        "Value" => StatsBase.median => "Value",
        "Value" => lower => "xmin",
        "Value" => upper => "xmax",
    )
    df."Affinity" = [hKav[hKav."IgG" .== r."Subclass", r."Cell"][1] for r in eachrow(df)]
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    pl1 = plot(
        df,
        x = "Affinity",
        y = "Value",
        ymin = "xmin",
        ymax = "xmax",
        color = "Cell",
        shape = "Subclass",
        Geom.point,
        Geom.errorbar,
        Scale.x_log10,
        Scale.y_log10,
        Guide.title("Recorded affinity vs. measured binding for single IgG"),
        Guide.xlabel("Recorded Affinity"),
        Guide.ylabel("Binding quantification"),
    )

    val_ratio = combine(groupby(df, ["Cell", "Subclass"])) do df
        (Ratio = df[df."Valency" .== Symbol("33"), "Value"][1] / df[df."Valency" .== Symbol("4"), "Value"][1],)
    end
    val_ratio."Affinity" = [hKav[hKav."IgG" .== r."Subclass", r."Cell"][1] for r in eachrow(val_ratio)]

    pl2 = plot(
        val_ratio,
        x = "Affinity",
        y = "Ratio",
        color = "Cell",
        shape = "Subclass",
        Geom.point,
        Scale.x_log10,
        Scale.y_log10,
        Guide.title("Recorded affinity vs. intervalency ratio"),
        Guide.xlabel("Recorded Affinity"),
        Guide.ylabel("33- to 4-valent binding quantification ratio"),
    )
    return pl1, pl2
end

igg_color_designation = Dict([humanIgG[i] => Scale.color_discrete().f(4)[i] for i = 1:length(humanIgG)])
igg_pair_color(iggA, iggB; tot = 5) = reverse([i for i in ColorScheme(range(igg_color_designation[iggA], igg_color_designation[iggB], length = tot))])

function plot_PCA_score(df; title = "Score")
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    layers = []
    for val in unique(df."Valency")
        for pair in unique(df."Subclass Pair")
            ddf = df[(df."Subclass Pair" .== pair) .& (df."Valency" .== val), :]
            sort!(ddf, ["%_2"])
            arrdf = DataFrame(xstart = Float64[], ystart = Float64[], xend = Float64[], yend = Float64[], Subclass = String[])
            for ii = 1:(nrow(ddf) - 1)
                push!(arrdf, [ddf[ii, "PC 1"], ddf[ii, "PC 2"], ddf[ii + 1, "PC 1"], ddf[ii + 1, "PC 2"], "Mixed"])
            end
            append!(layers, layer(arrdf, x = :xstart, y = :ystart, xend = :xend, yend = :yend, color = [colorant"black"], Geom.segment))
            # color=igg_pair_color(ddf."subclass_1"[1], ddf."subclass_2"[1]; tot=nrow(arrdf))
        end
    end

    df."Subclass" = copy(df."subclass_1")
    df[df."%_2" .== 1.0, "Subclass"] .= df[df."%_2" .== 1.0, "subclass_2"]
    df[(df."%_1" .< 1.0) .& (df."%_2" .< 1.0), "Subclass"] .= "Mixed"
    sdf = df[df."Subclass" .!= "Mixed", :]
    append!(layers, layer(df, x = "PC 1", y = "PC 2", color = [colorant"black"], size = [1mm], Geom.point))
    return plot(
        sdf,
        layers...,
        x = "PC 1",
        y = "PC 2",
        color = "Subclass",
        size = [3mm],
        Geom.point,
        Guide.title(title),
        Guide.xticks(ticks = [-2e4, -1e4, 0, 1e4], orientation = :horizontal),
        Guide.yticks(ticks = [-5e3, 0, 5e3]),
    )
end

function figure1()
    setGadflyTheme()

    p1, p2 = bindVSaff()

    # Specific IgG pair - receptor interaction
    df = averageMixData()
    igg12_1 = splot_origData(df[(df."Cell" .== "FcgRI") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG2"), :]; match_y = false)
    igg14_1 = splot_origData(df[(df."Cell" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG4"), :]; match_y = false)

    score, loading, vars_expl = mixtureDataPCA()
    vars = plot(
        DataFrame(Components = 1:length(vars_expl), R2X = vars_expl),
        x = "Components",
        y = "R2X",
        label = [@sprintf("%.2f %%", i * 100) for i in vars_expl],
        Geom.point,
        Geom.line,
        Geom.label,
        Scale.x_discrete,
        Scale.y_continuous(minvalue = 0.5),
        Guide.title("Variance Explained by PCA"),
        Guide.xlabel("Number of components"),
        Guide.ylabel("R2X"),
    )

    score_plot4 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs")
    score_plot33 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs")
    loading_plot = plot(loading, x = "PC 1", y = "PC 2", color = "Cell", label = "Cell", Geom.point, Geom.label, Guide.title("PCA Loadings"))


    pl = plotGrid((3, 4), [nothing, nothing, p1, p2, vars, score_plot4, score_plot33, loading_plot, igg12_1, igg14_1];)
    return draw(SVG("figure1.svg", 18inch, 12inch), pl)
end
