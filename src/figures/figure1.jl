""" Figure 1: show the mixture IC binding data """

using Printf

""" Original measurements with middle 50% as error bar """
function splot_origData(df)
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
        Scale.y_continuous(; maxvalue = ymax[cell]),
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

function plot_PCA_score(df)
    layers = []
    for val in unique(df."Valency")
        for pair in unique(df."Subclass Pair")
            append!(layers, layer(df[(df."Subclass Pair" .== pair) .& (df."Valency" .== val), :], x = "PC 1", y = "PC 2", color = [pair], Geom.line))
        end
    end
    return plot(df, layers..., x = "PC 1", y = "PC 2", color = "Subclass Pair", shape = "Valency", Geom.point, Guide.title("Score"))
end

function figure1()
    setGadflyTheme()

    df = averageMixData()
    pl1 = splot_origData(df[(df."Cell" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG3") .& (df."subclass_2" .== "IgG4"), :])
    pl2 = splot_origData(df[(df."Cell" .== "FcgRI") .& (df."subclass_1" .== "IgG3") .& (df."subclass_2" .== "IgG4"), :])
    pl3 = splot_origData(df[(df."Cell" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG4"), :])

    score_df, loading_df, vars_expl = mixtureDataPCA()
    vars = plot(
        x = 1:length(vars_expl),
        y = vars_expl,
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

    score_df[!, "Valency"] .= Symbol.(score_df[!, "Valency"])
    score_df."Subclass Pair" = score_df."subclass_1" .* "-" .* score_df."subclass_2"

    score_plot = plot_PCA_score(score_df)
    loading_plot = plot(loading_df, x = "PC 1", y = "PC 2", color = "Cell", label = "Cell", Geom.point, Geom.label, Guide.title("Loading"))

    pl = plotGrid((2, 4), [nothing, pl1, pl2, pl3, nothing, vars, score_plot, loading_plot]; widths = [1 1 1 1; 1 0.8 1.1 1.1])
    return draw(SVG("figure1.svg", 18inch, 8inch), pl)
end
