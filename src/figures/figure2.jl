""" Figure 2: we can accurately account for mixed ICs """

""" A general plotting function for Adjusted vs. Predicted plots """
function plotPredvsMeasured(
    df;
    xx = "Adjusted",
    yy = "Predict",
    xxlabel = "Actual",
    yylabel = "Predicted",
    color = "Cell",
    shape = "Valency",
    title = "Predicted vs Actual",
    R2pos = (4, 1),
)
    setGadflyTheme()

    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    df[(df[!, xx]) .< 1.0, xx] .= 1.0
    df[(df[!, yy]) .< 1.0, yy] .= 1.0

    r2 = R2((df[!, xx]), (df[!, yy]))
    return plot(
        df,
        x = xx,
        y = yy,
        xmin = "xmin",
        xmax = "xmax",
        color = color,
        shape = shape,
        Geom.point,
        "xmin" in names(df) ? Geom.errorbar : Guide.xlabel(xxlabel),
        Guide.ylabel(yylabel, orientation = :vertical),
        Guide.title(title),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(
            Scale.color_discrete().f(10)[1],
            Scale.color_discrete().f(10)[3],
            Scale.color_discrete().f(10)[2],
            Scale.color_discrete().f(10)[4:end]...,
        ),
        Geom.abline(color = "black"),
        Guide.annotation(compose(context(), text(R2pos[1], R2pos[2], "R<sup>2</sup> = " * @sprintf("%.4f", r2)), font("Helvetica-Bold"))),
        style(errorbar_cap_length = 0px),
    )
end

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

    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]

    pl = plot(
        layer(x = x, y = preds4, Geom.line, Theme(default_color = palette[1], line_width = 2px)),
        layer(x = x, y = preds33, Geom.line, Theme(default_color = palette[2], line_width = 2px)),
        Scale.x_continuous(labels = n -> "$IgGXname $(n*100)%\n$IgGYname $(100-n*100)%"),
        Guide.xlabel(""),
        Guide.ylabel("RFU", orientation = :vertical),
        Guide.xticks(orientation = :horizontal),
        Guide.title("$IgGXname-$IgGYname in $(df[1, "Cell"])"),
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

function figure2()
    data = loadMixData(; discard_small = true)
    res, df = fitMixMaster(data, fitKav = true)

    all_fit = plotPredvsMeasured(df; xx = "Value")
    pure_fit = plotPredvsMeasured(df[(df."%_1" .== 1) .| (df."%_2" .== 1), :]; xx = "Value", title = "Predicted vs Actual, single isotype only")

    p1 = splot_pred("FcgRIIIA-158F"; Lbound = true)
    p2 = splot_pred("FcgRIIIA-158F"; Lbound = false)

    draw(SVG("figure2.svg", 16inch, 9inch), plotGrid((2, 3), [nothing, all_fit, pure_fit, nothing, p1, p2]; sublabels = [1 1 1 0 1 1]))
end
