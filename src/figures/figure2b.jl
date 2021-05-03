""" Figure 2: we can accurately account for mixed ICs """

function figure2b()
    setGadflyTheme()
    Cellfit = false
    adjusted = true
    IgGx_Only = true

    data = loadMixData()

    if IgGx_Only
        data = data[data[!, "%_1"] .!= 10 / 100, :]
        data = data[data[!, "%_1"] .!= 33 / 100, :]
        data = data[data[!, "%_1"] .!= 66 / 100, :]
        data = data[data[!, "%_1"] .!= 90 / 100, :]
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

    draw(SVG("figure2b.svg", 700px, 600px), plotPredvsMeasured(df, "Value", r2))
end

function figure2c()
    pl = plotPredvsMeasured(PCA_dimred(), "PCA", "Predict", "Actual on imputed PC1")
    draw(SVG("figure2c.svg", 700px, 600px), pl)
end
