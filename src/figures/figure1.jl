""" Figure 1: show the mixture IC binding data """

function figure1()
    setGadflyTheme()
    draw(SVG("figure1.svg", 20inch, 16inch), plotMixOriginalData(loadMixData()))
    #draw(SVG("figure1A.svg", 20inch, 16inch), plotMixOriginalData(PCAData(; cutoff = 0.95)))

    df = loadMixData()
    df2 = predictMix(df)
    res = MixtureFit(df2; logscale = true)
    pl_pred = makeMixturePairSubPlots(res["df"])
    draw(SVG("figure1_all.svg", 20inch, 16inch), pl_pred)

    df3 = mixNormalExpBatch()
    df3 = df3[!, Not(["ymin", "ymax"])]
    df3[!, "Experiment"] .= "median"
    res = MixtureFit(df3; logscale = true)
    pl_pred = makeMixturePairSubPlots(res["df"])
    draw(SVG("figure1_median.svg", 20inch, 16inch), pl_pred)
end


""
