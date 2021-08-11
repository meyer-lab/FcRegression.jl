""" Figure 2: we can accurately account for mixed ICs """

""" A general plotting function for Adjusted vs. Predicted plots """
function plotPredvsMeasured(df; xx = "Adjusted", yy = "Predict", xxlabel = "Actual", 
    yylabel = "Predicted", color = "Valency", shape = "Cell")
    setGadflyTheme()

    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    df[(df[!, xx]) .< 1.0, xx] .= 1.0
    df[(df[!, yy]) .< 1.0, yy] .= 1.0

    xmins = "StdDev" in names(df) ? (df[!, xx] .- df[!, "StdDev"]) : df[!, xx]
    xmaxs = "StdDev" in names(df) ? (df[!, xx] .+ df[!, "StdDev"]) : xmins
    xmins[xmins .< 1.0] .= 1.0
    xmaxs[xmaxs .< 1.0] .= 1.0

    r2 = R2((df[!, xx]), (df[!, yy]))
    return plot(
        df,
        x = xx,
        y = yy,
        xmin = xmins,
        xmax = xmaxs,
        color = color,
        shape = shape,
        Geom.point,
        "StdDev" in names(df) ? Geom.errorbar : Guide.xlabel(xxlabel),
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

    if !("Adjusted" in names(data))
        data[!, "Adjusted"] .= data[!, "Value"]
    end
    df = predictMix(data)

    draw(SVG("figure2.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df; xx = "Value")]))
end

