function figure4()

    setGadflyTheme()
    
    dataType = "ITP"
    df = importDepletion(dataType)
    color = (dataType == "HIV") ? "Label" : "Background"
    shape = "Condition"
    L0 = 1e-9
    f = 4
    murine = true

    res, odf, Cell_df, ActI_df = regressionResult(dataType; L0 = L0, f = f, murine = murine)
    @assert all(in(names(odf)).([color, shape]))

    ITP_Dep = plotDepletionSynergy(IgGX = 2, IgGY = 4; L0 = L0, f = f, murine = murine, c1q = ("C1q" in Cell_df.Component), dataType = "ITP", fit = res, neutralization = ("Neutralization" in names(df)))
    ITP_Dep_All = plotSynergy(L0, f; murine = murine, dataType = dataType, fit = res)
    ITP_Kupffer_Act = plotDepletionSynergy(IgGX = 2, IgGY = 4; L0 = L0, f = f, murine = murine, c1q = ("C1q" in Cell_df.Component), dataType = "ITP", fit = res, Cellidx = 6, neutralization = ("Neutralization" in names(df)))
    ITP_Kupffer_Act_All = plotSynergy(L0, f; murine = murine, dataType = dataType, fit = res, Cellidx = 6)
    ITP_Kupffer_Bound = plotDepletionSynergy(IgGX = 2, IgGY = 4; L0 = L0, f = f, murine = murine, c1q = ("C1q" in Cell_df.Component), dataType = "ITP", fit = res, Cellidx = 6, Rbound = true, neutralization = ("Neutralization" in names(df)))
    ITP_Kupffer_Bound_All = plotSynergy(L0, f; murine = murine, dataType = dataType, fit = res, Cellidx = 6, Rbound = true)

    draw(
        SVG("figure4.svg", 10inch, 8inch),
        plotGrid((2, 3), [ITP_Dep ITP_Kupffer_Act ITP_Kupffer_Bound ITP_Dep_All ITP_Kupffer_Act_All ITP_Kupffer_Bound_All]; widths = [3 3 3; 4 4 4]),
    )
end
