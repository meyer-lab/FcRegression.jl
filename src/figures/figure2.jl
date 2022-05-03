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
    R2pos = (3, 1),
    clip2one = true,
)
    setGadflyTheme()

    df = deepcopy(df)
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    if clip2one
        df[(df[!, xx]) .< 1.0, xx] .= 1.0
        df[(df[!, yy]) .< 1.0, yy] .= 1.0
    end

    r2 = R2((df[!, xx]), (df[!, yy]))
    println(r2)
    errbar = "xmin" in names(df)
    return plot(
        df,
        x = xx,
        y = yy,
        xmin = (errbar ? "xmin" : xx),
        xmax = (errbar ? "xmax" : xx),
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
        Guide.annotation(
            compose(
                context(),
                text(R2pos[1], R2pos[2], "<i>R</i><sup>2</sup> = " * @sprintf("%.4f", r2)),
                stroke("black"),
                fill("black"),
                font("Helvetica-Bold"),
            ),
        ),
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

function figure2()
    data = loadMixData()
    data = averageMixData(data)
    raw_predict = predictMix(data)

    raw_pred_pl = plotPredvsMeasured(raw_predict; xx = "Value", xxlabel = "Measured", title = "Prediction without fitting", R2pos = (3, 1))


    c = runMCMC("MCMC_nuts_0502.dat")
    df = loadMixData()
    pl1 = MCMC_params_predict_plot(c, df; xx = "Value", yy = "Predict", title = "All predictions with \nsingle IgG fitted params")
    pl2 = MCMC_params_predict_plot(
        c,
        df[(df."%_1" .!= 1.0) .& (df."%_2" .!= 1.0), :];
        xx = "Value",
        yy = "Predict",
        title = "Mixture predictions with \nsingle IgG fitted params",
    )
    _, pl_igg, _ = plot_MCMC_affinity(c)
    rob1, rob2 = validateRobinett("MCMC_robinett_0503.dat", c; mcmc_iter = 100)

    pp = plotGrid((4, 3), [nothing, raw_pred_pl, pl1, pl2, pl_igg[1], pl_igg[2], pl_igg[3], pl_igg[4], rob1, rob2])
    draw(PDF("figure2.pdf", 10inch, 12inch), pp)
end
