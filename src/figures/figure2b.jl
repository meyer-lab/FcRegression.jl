""" Figure 2: we can accurately account for mixed ICs """

function figure2b()
    data = loadMixData()
    data = data[data[!, "%_1"] .!= 10 / 100, :]
    data = data[data[!, "%_1"] .!= 33 / 100, :]
    data = data[data[!, "%_1"] .!= 66 / 100, :]
    data = data[data[!, "%_1"] .!= 90 / 100, :]

    df = MixtureCellSeparateFit(data; logscale = true, adjusted = false)
    draw(SVG("figure2b.svg", 700px, 600px), plotPredvsMeasured(df, "Value"))
end

function figure2c()
    pl = plotPredvsMeasured(PCA_dimred(), "PCA", "Predict", "Actual on imputed PC1")
    draw(SVG("figure2c.svg", 700px, 600px), pl)
end
