""" Figure 2: we can accurately account for mixed ICs """

function figure2b()
    setGadflyTheme()

    data = loadMixData()
    data = data[data[!,"%_1"].!=0,:]
    data = data[data[!,"%_1"].!=10,:]
    data = data[data[!,"%_1"].!=33,:]
    data = data[data[!,"%_1"].!=66,:]
    data = data[data[!,"%_1"].!=90,:]
    data[!, "%_1"] ./= 100.0
    df = MixtureCellSeparateFit(data; logscale = true, adjusted = false)
    df[!, "Value"] = abs.(df[!, "Value"])
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    pl = plot(
        df,
        x = "Value",
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

    draw(SVG("figure2b.svg", 1300px, 600px), plotGrid((1, 2), [nothing, pl]))
end
