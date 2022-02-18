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
)
    setGadflyTheme()

    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    df[(df[!, xx]) .< 1.0, xx] .= 1.0
    df[(df[!, yy]) .< 1.0, yy] .= 1.0

    r2 = R2((df[!, xx]), (df[!, yy]))
    R2anno = "R<sup>2</sup> = " * @sprintf("%.4f", r2)
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
        Guide.annotation(compose(context(), text(160, 60000, R2anno), fill("black"), fontsize(10pt), font("Helvetica"))),
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
    raw_predict = mixturePredictions(averageMixData(data))

    raw_pred_pl = plotPredvsMeasured(raw_predict; xx = "Value", xxlabel = "Measured", title = "Prediction without fitting")

    reg_res, reg_fitdf = fitMixMaster(data; fitKav = false)
    reg_allPL = plotPredvsMeasured(reg_fitdf; xx = "Value", title = "Fit all except Kav, all")
    reg_onePL = plotPredvsMeasured(
        reg_fitdf[(reg_fitdf."%_1" .== 1) .| (reg_fitdf."%_2" .== 1), :];
        xx = "Value",
        title = "Fit all except Kav, single isotypes",
    )

    res = exp.(reg_res.minimizer)
    _, kfitdf = fitMixMaster(data; fitRVX = false, recepExp = res[1:(end - 3)], vals = res[(end - 2):(end - 1)], KxStar = res[end], fitKav = true)
    kfit_allPL = plotPredvsMeasured(kfitdf; xx = "Value", title = "Fit Kav, all")
    kfit_onePL = plotPredvsMeasured(kfitdf[(kfitdf."%_1" .== 1) .| (kfitdf."%_2" .== 1), :]; xx = "Value", title = "Fit Kav, single isotypes")

    p1 = splot_pred("FcgRIIIA-158F"; Lbound = true)
    p2 = splot_pred("FcgRIIIA-158F"; Lbound = false)

    draw(
        SVG("figure2.svg", 16inch, 13inch),
        plotGrid((3, 3), [nothing, raw_pred_pl, reg_allPL, reg_onePL, kfit_allPL, kfit_onePL, p1, p2, nothing]; sublabels = [1 1 1 1 1 1 1 1 0]),
    )
end
