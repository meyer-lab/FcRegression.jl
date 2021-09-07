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

""" Individual measurement with prediction curve """
function splot_contPred(df; logscale = false)
    df = copy(df)
    @assert length(unique(df."Cell")) == 1
    @assert length(unique(df."subclass_1")) == 1
    @assert length(unique(df."subclass_2")) == 1
    IgGXname = unique(df."subclass_1")[1]
    IgGYname = unique(df."subclass_2")[1]

    x = 0:0.01:1
    df4 = df[(df."Valency" .== 4), :]
    preds4 = [predictMix(df4[1, :], IgGXname, IgGYname, i, 1 - i) for i in x]
    df33 = df[(df."Valency" .== 33), :]
    preds33 = [predictMix(df33[1, :], IgGXname, IgGYname, i, 1 - i) for i in x]

    if !("Adjusted" in names(df))
        df[!, "Adjusted"] .= df[!, "Value"]
    end
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]
    pl = plot(
        layer(x = x, y = preds4, Geom.line, Theme(default_color = palette[1], line_width = 2px)),
        layer(x = x, y = preds33, Geom.line, Theme(default_color = palette[2], line_width = 2px)),
        layer(df, x = "%_1", y = "Adjusted", color = "Valency", shape = "Experiment"),
        Scale.x_continuous(labels = n -> "$IgGXname $(n*100)%\n$IgGYname $(100-n*100)%"),
        (logscale ? Scale.y_log10(minvalue = 1, maxvalue = 1e6) : Scale.y_continuous),
        Scale.color_discrete_manual(palette[1], palette[2]),
        Guide.xlabel(""),
        Guide.ylabel("RFU", orientation = :vertical),
        Guide.xticks(orientation = :horizontal),
        Guide.title("$IgGXname-$IgGYname in $(df[1, "Cell"])"),
    )
    return pl
end

function plot_PCA_score(df; title="Score")
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    layers = []
    for val in unique(df."Valency")
        for pair in unique(df."Subclass Pair")
            append!(layers, layer(df[(df."Subclass Pair" .== pair) .& (df."Valency" .== val), :], x = "PC 1", y = "PC 2", color = [pair], Geom.line))
        end
    end
    return plot(df, layers..., x = "PC 1", y = "PC 2", color = "Subclass Pair", Geom.point, Guide.title(title))
end

function figure1()
    setGadflyTheme()

    df = averageMixData()
    pl1 = splot_origData(df[(df."Cell" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG3") .& (df."subclass_2" .== "IgG4"), :]; match_y = false)
    pl2 = splot_origData(df[(df."Cell" .== "FcgRI") .& (df."subclass_1" .== "IgG2") .& (df."subclass_2" .== "IgG4"), :]; match_y = false)
    pl3 = splot_origData(df[(df."Cell" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG4"), :]; match_y = false)

    _, _, vars_expl = mixtureDataPCA()
    score_df4, loading_df4, vars_expl4 = mixtureDataPCA(; val = 4)
    score_df33, loading_df33, vars_expl33 = mixtureDataPCA(; val = 33)
    R2Xdat = vcat([DataFrame(Components = 1:length(vars_expl), Data = i, R2X = j) 
        for (i, j) in [("Overall", vars_expl), ("4-valent IC", vars_expl4), ("33-valent IC", vars_expl33)]]...)
    vars = plot(
        R2Xdat,
        x = "Components",
        y = "R2X",
        color = "Data",
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

    score_plot4 = plot_PCA_score(score_df4; title="Score, 4-valent ICs")
    loading_plot4 = plot(loading_df4, x = "PC 1", y = "PC 2", color = "Cell", label = "Cell", Geom.point, Geom.label, Guide.title("Loading, 4-valent ICs"))
    score_plot33 = plot_PCA_score(score_df33; title="Score, 33-valent ICs")
    loading_plot33 = plot(loading_df33, x = "PC 1", y = "PC 2", color = "Cell", label = "Cell", Geom.point, Geom.label, Guide.title("Loading, 33-valent ICs"))

    pl = plotGrid(
        (3, 4),
        [nothing, pl1, pl2, pl3, nothing, vars, score_plot4, loading_plot4, score_plot33, loading_plot33];
        sublabels = [1 1 1 1 0 1 1 1 1 1 0 0],
    )
    return draw(SVG("figure1.svg", 18inch, 12inch), pl)
end
