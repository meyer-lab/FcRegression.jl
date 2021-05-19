""" Figure 2: we can accurately account for mixed ICs """

function figure2b()
    Cellfit = true
    adjusted = true
    IgGx_Only = true

    data = loadMixData(avg=true)

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

    draw(SVG("figure2b.svg", 700px, 600px), plotPredvsMeasured(df; xx = xvar))
end

function figure2c()
    pl = plotPredvsMeasured(PCA_dimred(), xx="PCA", yy="Predict", xxlabel="Actual on imputed PC1")
    draw(SVG("figure2c.svg", 700px, 600px), pl)
end
