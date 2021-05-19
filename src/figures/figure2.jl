""" Figure 2: we can accurately account for mixed ICs """

function plotPredvsMeasured(df; xx = "Adjusted", yy = "Predict", 
        xxlabel = "Actual", yylabel = "Predicted", color = "Valency", shape = "Cell")
    setGadflyTheme()
    df[!, color] .= Symbol.(df[!, color])
    df[!, shape] .= Symbol.(df[!, shape])
    df[(df[!, xx]) .< 1.0, xx] .= 1.0
    df[(df[!, yy]) .< 1.0, yy] .= 1.0
    

    r2 = R2((df[!, xx]), (df[!, yy]))

    return plot(
        df,
        x = xx,
        y = yy,
        color = color,
        shape = shape,
        Geom.point,
        Guide.xlabel(xxlabel),
        Guide.ylabel(yylabel, orientation = :vertical),
        Guide.title("R^2: $r2"),
        Scale.x_log10,
        Scale.y_log10,
        Scale.color_discrete_manual(Scale.color_discrete().f(10)[1], Scale.color_discrete().f(10)[3], 
            Scale.color_discrete().f(10)[2], Scale.color_discrete().f(10)[4:end]...),
        Geom.abline(color = "green"),
    )
end


function figure2()
    df = MixtureCellSeparateFit(loadMixData(avg=true); logscale = true)
    draw(SVG("figure2.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df)]))
end
