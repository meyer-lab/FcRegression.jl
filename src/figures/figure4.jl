""" Figure 4: validating the in vivo effect of mixtures """

function figure4()

    setGadflyTheme()

    dataType = "ITP"
    df = importDepletion(dataType)
    color = "Label"
    shape = "Condition"
    L0 = 1e-9
    f = 4

    res, odf, Cell_df, ActI_df = regressionResult(dataType; L0 = L0, f = f)

    ITP_Dep = plotDepletionSynergy(1, 2; L0 = L0, f = f, dataType = "ITP", fit = res)
    ITP_Dep_All = plotSynergy(L0, f; dataType = dataType, fit = res)
    ITP_Kupffer_Act = plotDepletionSynergy(1, 2; L0 = L0, f = f, dataType = "ITP", fit = res, Cellidx = 6)
    ITP_Kupffer_Act_All = plotSynergy(L0, f; dataType = dataType, fit = res, Cellidx = 6)
    ITP_Kupffer_Bound = plotDepletionSynergy(1, 2; L0 = L0, f = f, dataType = "ITP", fit = res, Cellidx = 6, Rbound = true)
    ITP_Kupffer_Bound_All = plotSynergy(L0, f; dataType = dataType, fit = res, Cellidx = 6, Rbound = true)

    draw(
        SVG("figure4.svg", 10inch, 8inch),
        plotGrid((2, 3), [ITP_Dep ITP_Kupffer_Act ITP_Kupffer_Bound ITP_Dep_All ITP_Kupffer_Act_All ITP_Kupffer_Bound_All]; widths = [3 3 3; 4 4 4]),
    )
end
