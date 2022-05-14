""" Figure 1: show the mixture IC binding data """

using ColorSchemes

""" Original measurements with middle 50% as error bar """
function splot_origData(df; match_y = true)
    cell = unique(df."Receptor")[1]
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
    df = df[!, ["Valency", "Receptor", "Subclass", "Value"]]
    df = combine(
        groupby(df, ["Valency", "Receptor", "Subclass"]),
        "Value" => StatsBase.median => "Value",
        "Value" => lower => "xmin",
        "Value" => upper => "xmax",
    )
    df."Affinity" = [hKav[hKav."IgG" .== r."Subclass", r."Receptor"][1] for r in eachrow(df)]
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    pl1 = plot(
        df,
        x = "Affinity",
        y = "Value",
        ymin = "xmin",
        ymax = "xmax",
        color = "Receptor",
        shape = "Subclass",
        Geom.point,
        Geom.errorbar,
        Scale.x_log10,
        Scale.y_log10,
        Guide.title("Recorded affinity vs. measured binding for single IgG"),
        Guide.xlabel("Recorded Affinity"),
        Guide.ylabel("Binding quantification"),
    )

    val_ratio = combine(groupby(df, ["Receptor", "Subclass"])) do df
        (Ratio = df[df."Valency" .== Symbol("33"), "Value"][1] / df[df."Valency" .== Symbol("4"), "Value"][1],)
    end
    val_ratio."Affinity" = [hKav[hKav."IgG" .== r."Subclass", r."Receptor"][1] for r in eachrow(val_ratio)]

    pl2 = plot(
        val_ratio,
        x = "Affinity",
        y = "Ratio",
        color = "Receptor",
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

function plot_PCA_score(df; title = "Score", xx = "PC 1", yy = "PC 2")
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    layers = []
    for val in unique(df."Valency")
        for pair in unique(df."Subclass Pair")
            ddf = df[(df."Subclass Pair" .== pair) .& (df."Valency" .== val), :]
            sort!(ddf, ["%_2"])
            arrdf = DataFrame(xstart = Float64[], ystart = Float64[], xend = Float64[], yend = Float64[], Subclass = String[])
            for ii = 1:(nrow(ddf) - 1)
                push!(arrdf, [ddf[ii, xx], ddf[ii, yy], ddf[ii + 1, xx], ddf[ii + 1, yy], "Mixed"])
            end
            append!(layers, layer(arrdf, x = :xstart, y = :ystart, xend = :xend, yend = :yend, color = [colorant"black"], Geom.segment))
            # color=igg_pair_color(ddf."subclass_1"[1], ddf."subclass_2"[1]; tot=nrow(arrdf))
        end
    end

    df."Subclass" = copy(df."subclass_1")
    df[df."%_2" .== 1.0, "Subclass"] .= df[df."%_2" .== 1.0, "subclass_2"]
    df[(df."%_1" .< 1.0) .& (df."%_2" .< 1.0), "Subclass"] .= "Mixed"
    sdf = df[df."Subclass" .!= "Mixed", :]
    append!(layers, layer(df, x = xx, y = yy, color = [colorant"black"], size = [1mm], Geom.point))
    return plot(
        sdf,
        layers...,
        x = xx,
        y = yy,
        color = "Subclass",
        size = [3mm],
        Geom.point,
        Guide.title(title),
        Guide.xticks(ticks = [-1.5e4, 0, 1.5e4], orientation = :horizontal),
        Guide.yticks(ticks = [-6e3, 0, 6e3]),
    )
end

function figure1()
    setGadflyTheme()

    p1, p2 = bindVSaff()

    # Specific IgG pair - receptor interaction
    df = averageMixData()
    igg12_1 = splot_origData(df[(df."Receptor" .== "FcgRI") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG2"), :]; match_y = false)
    igg14_1 = splot_origData(df[(df."Receptor" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG4"), :]; match_y = false)

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

    SP4 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs", xx = "PC 1", yy = "PC 2")
    SP33 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs", xx = "PC 1", yy = "PC 2")
    SP4_13 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs", xx = "PC 1", yy = "PC 3")
    SP33_13 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs", xx = "PC 1", yy = "PC 3")
    LP = plot(
        loading,
        x = "PC 1",
        y = "PC 2",
        color = "Receptor",
        label = "Receptor",
        Geom.point,
        Geom.label,
        Guide.title("PCA Loadings"),
        Scale.x_continuous(minvalue = -1.0, maxvalue = 1.0),
        Scale.y_continuous(minvalue = -1.0, maxvalue = 1.0),
    )
    LP_13 = plot(
        loading,
        x = "PC 1",
        y = "PC 3",
        color = "Receptor",
        label = "Receptor",
        Geom.point,
        Geom.label,
        Guide.title("PCA Loadings"),
        Scale.x_continuous(minvalue = -1.0, maxvalue = 1.0),
        Scale.y_continuous(minvalue = -1.0, maxvalue = 1.0),
    )

    pl = plotGrid((3, 5), [nothing, nothing, p1, p2, vars, SP4, SP33, LP, igg12_1, igg14_1, SP4_13, SP33_13, LP_13];)
    return draw(PDF("figure1.pdf", 20inch, 10inch), pl)
end
