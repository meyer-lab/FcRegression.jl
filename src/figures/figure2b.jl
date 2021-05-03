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
        xvar = "Adjusted"
        df[(df[!, "Adjusted"]) .< 1.0, "Adjusted"] .= 1.0
    elseif adjusted
        df = MixtureFit(data; logscale = true)["df"]
        xvar = "Adjusted"
        df[(df[!, "Adjusted"]) .< 1.0, "Adjusted"] .= 1.0
    else
        df = predictMix(data)
        xvar = "Value"
        df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0
    end

    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    r2 = R2((df[!, xvar]), (df[!, "Predict"]))

    draw(SVG("figure2b.svg", 700px, 600px), plotPredvsMeasured(df, "Value", r2))
end

function figure2c()
    pl = plotPredvsMeasured(PCA_dimred(), "PCA", "Predict", "Actual on imputed PC1")
    draw(SVG("figure2c.svg", 700px, 600px), pl)
end
