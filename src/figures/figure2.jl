""" Figure 2: we can accurately account for mixed ICs """

function figure2()
    setGadflyTheme()

    df = MixtureCellSeparateFit(loadMixData(); logscale = true)
    df[!, "%_1"] ./= 100.0
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    pl = plot(
        df,
        x = "Adjusted",
        y = "Predict",
        color = "Valency",
        shape = "Cell",
        Geom.point,
        Guide.xlabel("Actual"),
        Guide.ylabel("Predicted", orientation = :vertical),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]),
        Geom.abline(color = "green"),
    )

    draw(SVG("figure2.svg", 1300px, 600px), plotGrid((1, 2), [nothing, pl]))
end
