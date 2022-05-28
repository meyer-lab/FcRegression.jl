""" Figure 1: show the mixture IC binding data """

## Disabled because current predMix() cannot handle
## TODO: enable predMix() to take this
#=
""" Individual measurement with prediction curve """
function splot_contPred(df)
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

    ymax =
        Dict("FcgRI" => 5e4, "FcgRIIA-131H" => 2e5, "FcgRIIA-131R" => 5e4, "FcgRIIB-232I" => 600, "FcgRIIIA-158F" => 1.5e5, "FcgRIIIA-158V" => 2.5e5)
    cell = df[1, "Cell"]
    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]

    pl = plot(
        layer(x = x, y = preds4, Geom.line, Theme(default_color = palette[1], line_width = 2px)),
        layer(x = x, y = preds33, Geom.line, Theme(default_color = palette[2], line_width = 2px)),
        Scale.x_continuous(labels = n -> "$IgGXname $(n*100)%\n$IgGYname $(100-n*100)%"),
        Scale.y_continuous(; minvalue = 0.0, maxvalue = ymax[cell]),
        Guide.xlabel(""),
        Guide.ylabel("RFU", orientation = :vertical),
        Guide.xticks(orientation = :horizontal),
        Guide.title("$IgGXname-$IgGYname in $cell"),
        Guide.manual_color_key("Valency", ["4", "33"], [palette[1], palette[2]]),
    )
    return pl
end

""" Prediction curve """
function splot_pred(cell; Lbound = true)
    x = 0:0.01:1
    preds = [[predictMix(cell, 33, IgGXname, "IgG2", i, 1 - i; Lbound = Lbound) for i in x] for IgGXname in ["IgG1", "IgG3", "IgG4"]]
    palette = Scale.color_discrete().f(10)
    pl = plot(
        layer(x = x, y = preds[1], Geom.line, Theme(default_color = palette[1], line_width = 2px)),
        layer(x = x, y = preds[2], Geom.line, Theme(default_color = palette[2], line_width = 2px)),
        layer(x = x, y = preds[3], Geom.line, Theme(default_color = palette[3], line_width = 2px)),
        Scale.x_continuous(labels = n -> "IgGX $(n*100)%\n IgG2 $(100-n*100)%"),
        Scale.y_log10(minvalue = 1),
        Guide.manual_color_key("Subclass Pair", ["IgG1-IgG2", "IgG3-IgG2", "IgG4-IgG2"], [palette[1], palette[2], palette[3]]),
        Guide.xlabel(""),
        Guide.ylabel("RFU", orientation = :vertical),
        Guide.xticks(orientation = :horizontal),
        Guide.title("Predicted " * (Lbound ? "binding" : "multimerization") * " to $cell"),
    )
    return pl
end
=#

function figureS1()
    setGadflyTheme()
    draw(PDF("figureS1.pdf", 20inch, 16inch), plotMixSubplots(splot_origData, averageMixData()))
    draw(PDF("figureS3.pdf", 20inch, 16inch), plotMixSubplots(splot_contPred, loadMixData()))
end
