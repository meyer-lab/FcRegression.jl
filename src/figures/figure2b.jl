""" Figure 2: we can accurately account for mixed ICs """

function figure2b()
    setGadflyTheme()

    data = loadMixData()
    data = data[data[!, "%_1"] .!= 0, :]
    data = data[data[!, "%_1"] .!= 10, :]
    data = data[data[!, "%_1"] .!= 33, :]
    data = data[data[!, "%_1"] .!= 66, :]
    data = data[data[!, "%_1"] .!= 90, :]
    data[!, "%_1"] ./= 100.0

    df = predictMix(data)

    #=adjusted = false
    Cellfit = true
    if Cellfit
        df = MixtureCellSeparateFit(data; logscale = true, adjusted = adjusted)
    else
        dict = MixtureFit(data; logscale = true, adjusted = adjusted)
        df = dict["df"]
    end
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    display(df)

    if adjusted
        xval = "Adjusted"
        df[(df[!, "Adjusted"]) .< 1.0, "Adjusted"] .= 1.0
        r2 = R2((df[!, "Adjusted"]), (df[!, "Predict"]))
    else
        xval = "Value"
        df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0
        r2 = R2((df[!, "Value"]), (df[!, "Predict"]))
    end=#

    pl = plot(
        df,
        x = "Value",
        y = "Predict",
        color = "Valency",
        shape = "Cell",
        Geom.point,
        Guide.xlabel("Actual"),
        Guide.ylabel("Predicted", orientation = :vertical),
        #Guide.title("R^2: $r2"),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]),
        Geom.abline(color = "green"),
    )

    draw(SVG("figure2b.svg", 1300px, 600px), plotGrid((1, 2), [nothing, pl]))
end
