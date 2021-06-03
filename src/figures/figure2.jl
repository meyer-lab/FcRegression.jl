""" Figure 2: we can accurately account for mixed ICs """

function plotPredvsMeasured(df; xx = "Adjusted", yy = "Predict", 
    xxlabel = "Actual", yylabel = "Predicted", color = "Valency", shape = "Cell")
    setGadflyTheme()

    df[!, "Valency"] .= Symbol.(df[!, "Valency"])
    df[(df[!, xx]) .< 1.0, xx] .= 1.0
    df[(df[!, yy]) .< 1.0, yy] .= 1.0

    xmins, xmaxs = errorBars(df; xx)

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
        "StdDev" in names(df) ? Geom.errorbar : ,
        Guide.xlabel(xxlabel),
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


function figure2(Cellfit = true, adjusted = true, IgGx_Only = false, avg = true)

    if Cellfit == true && adjusted == false
        @assert (Cellfit === adjusted) "Adjusted must be true if Cellfit is true"
    end

    if avg
        data = averageData(loadMixData())
    else
        data = loadMixData()
    end

    if IgGx_Only
        data = data[data[!, "%_1"] .!= 10 / 100, :]
        data = data[data[!, "%_1"] .!= 33 / 100, :]
        data = data[data[!, "%_1"] .!= 66 / 100, :]
        data = data[data[!, "%_1"] .!= 90 / 100, :]
    end

    if Cellfit
        df = MixtureCellSeparateFit(data; logscale = true)
        xvar = "Adjusted"
    elseif adjusted
        df = MixtureFit(data; logscale = true)["df"]
        xvar = "Adjusted"
    else
        df = predictMix(data)
        xvar = "Value"
    end

    draw(SVG("figure2.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df; xx = xvar)]))
end

function figure2c()
    pl = plotPredvsMeasured(PCA_dimred(); xx="PCA", yy="Predict", xxlabel="Actual on imputed PC1")
    draw(SVG("figure2c.svg", 700px, 600px), pl)
end
