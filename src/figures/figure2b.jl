""" Figure 2: we can accurately account for mixed ICs """

function figure2b()
    setGadflyTheme()
    Cellfit = false
    adjusted = false
    IgGx_Only = true

    data = loadMixData()

    if IgGx_Only
        data[data[!, "%_1"] .<= 1.0, "%_1"] .= 0
        data = data[data[!, "%_1"] .!= 0, :]
    end

    if Cellfit
        df = MixtureCellSeparateFit(data; logscale = true)
    elseif adjusted
        dict = MixtureFit(data; logscale = true)
        df = dict["df"]
    else
        df = predictMix(data)
    end

    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    if adjusted || Cellfit
        xvar = "Adjusted"
        df[(df[!, "Adjusted"]) .< 1.0, "Adjusted"] .= 1.0
    else
        xvar = "Value"
        df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0
    end

    r2 = R2((df[!, xvar]), (df[!, "Predict"]))
    display(df)

    pl = plot(
        df,
        x = xvar,
        y = "Predict",
        color = "Valency",
        shape = "Cell",
        Geom.point,
        Guide.xlabel("Actual"),
        Guide.ylabel("Predicted", orientation = :vertical),
        Guide.title("R^2: $r2"),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]),
        Geom.abline(color = "green"),
    )

    draw(SVG("figure2b.svg", 1300px, 600px), plotGrid((1, 2), [nothing, pl]))
end
