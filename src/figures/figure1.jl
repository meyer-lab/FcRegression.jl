""" Figure 1: show the mixture IC binding data """

using Printf

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
        Guide.title("Affinity vs Binding"),
        Guide.xlabel("Affinity"),
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
        Guide.title("Affinity vs Ratio"),
        Guide.xlabel("Affinity"),
        Guide.ylabel("Valency 33 to 4 binding quantification ratio"),
    )
    return pl1, pl2
end

function plot_PCA_score(df; title = "Score")
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    layers = []
    for val in unique(df."Valency")
        for pair in unique(df."Subclass Pair")
            append!(layers, layer(df[(df."Subclass Pair" .== pair) .& (df."Valency" .== val), :], x = "PC 1", y = "PC 2", color = [pair], Geom.line))
        end
    end
    return plot(
        df,
        layers...,
        x = "PC 1",
        y = "PC 2",
        color = "Subclass Pair",
        Geom.point,
        Guide.title(title),
        Guide.xticks(ticks = [-2e4, -1e4, 0, 1e4], orientation = :horizontal),
        Guide.yticks(ticks = [-5e3, 0, 5e3]),
    )
end

function figure1()
    setGadflyTheme()

    p1, p2 = bindVSaff()

    df = averageMixData()
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

    score_plot4 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "Score, 4-valent ICs")
    score_plot33 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "Score, 33-valent ICs")
    loading_plot = plot(loading, x = "PC 1", y = "PC 2", color = "Cell", label = "Cell", Geom.point, Geom.label, Guide.title("Loading"))

    # Specific IgG pair - receptor interaction
    igg12_1 = splot_origData(df[(df."Cell" .== "FcgRI") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG2"), :]; match_y = false)
    igg14_1 = splot_origData(df[(df."Cell" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG4"), :]; match_y = false)
    pl = plotGrid(
        (3, 4),
        [nothing, p1, p2, nothing, nothing, igg12_1, igg14_1, nothing, vars, score_plot4, score_plot33, loading_plot];
        sublabels = [1 1 1 1 0 1 1 0 1 1 1 1],
    )
    return draw(SVG("figure1.svg", 18inch, 12inch), pl)
end
