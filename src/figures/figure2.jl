""" Figure 2: we can accurately account for mixed ICs """

function plotPredvsMeasured(df, xx = "Adjusted", yy = "Predict", xxlabel = "Actual", yylabel = "Predicted")
    setGadflyTheme()
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    df[(df[!, xx]) .< 1.0, xx] .= 1.0
    df[(df[!, yy]) .< 1.0, yy] .= 1.0

    return plot(
        df,
        x = xx,
        y = yy,
        color = "Valency",
        shape = "Cell",
        Geom.point,
        Guide.xlabel(xxlabel),
        Guide.ylabel(yylabel, orientation = :vertical),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]),
        Geom.abline(color = "green"),
    )
end


function figure2()
    df = MixtureCellSeparateFit(loadMixData(); logscale = true)
    draw(SVG("figure2.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df)]))
end
