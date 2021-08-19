""" Figure 2: we can accurately account for mixed ICs """

""" A general plotting function for Adjusted vs. Predicted plots """
function plotPredvsMeasured(df; xx = "Adjusted", yy = "Predict", xxlabel = "Actual", yylabel = "Predicted", color = "Valency", shape = "Cell")
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
        Guide.title("R^2: $r2"),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(
            Scale.color_discrete().f(10)[1],
            Scale.color_discrete().f(10)[3],
            Scale.color_discrete().f(10)[2],
            Scale.color_discrete().f(10)[4:end]...,
        ),
        Geom.abline(color = "green"),
    )
end


function figure2(IgGx_Only = false)
    data = loadMixData()

    if IgGx_Only  # only one IgG subclass
        data = data[(data[!, "%_1"] .== 1.0) .| (data[!, "%_1"] .== 0.0), :]
    end
    df = fitMixMaster(data)[2]

    draw(SVG("figure2.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df; xx = "Value")]))
end
