""" Figure 5: analyze synergism in IgG mixture """

function figure5()

    setGadflyTheme()

    dataType = "ITP"
    df = importDepletion(dataType)
    color = "Label"
    shape = "Condition"
    L0 = 1e-9
    f = 4

    res, _, _ = regressionResult(dataType; L0 = L0, f = f, murine = true)

    Xdf = modelPred(df; L0 = L0, f = f, murine = true, cellTypes = nothing)
    res = fitRegNNLS(Xdf; murine = murine, cellTypes = nothing, exp_method = exp_method)

    ITP_Kupffer_Bound = plotDepletionSynergy(1, 2; L0 = L0, f = f, dataType = "ITP", fit = res, Cellidx = 6, Recepidx = 2, Rbound = true)
    ITP_Kupffer_Bound_All = plotSynergy(L0, f; murine = true, dataType = dataType, fit = res, Cellidx = 6, Recepidx = 2, Rbound = true)
    ITP_Kupffer_Act = plotDepletionSynergy(1, 2; L0 = L0, f = f, dataType = "ITP", fit = res, Cellidx = 6)
    ITP_Kupffer_Act_All = plotSynergy(L0, f; murine = true, dataType = dataType, fit = res, Cellidx = 6)
    ITP_Dep = plotDepletionSynergy(1, 2; L0 = L0, f = f, dataType = "ITP", fit = res)
    ITP_Dep_All = plotSynergy(L0, f; murine = true, dataType = dataType, fit = res)


    draw(
        SVG("figure5.svg", 10inch, 8inch),
        plotGrid((2, 3), [ITP_Kupffer_Bound ITP_Kupffer_Act ITP_Dep ITP_Kupffer_Bound_All ITP_Kupffer_Act_All ITP_Dep_All]),
    )
end
