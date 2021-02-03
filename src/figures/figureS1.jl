""" Figure S1 """

function figureS1()

    setGadflyTheme()

    dataType = "ITP"
    df = importDepletion(dataType)
    L0 = 1e-9
    f = 4

    res, odf, Cell_df, ActI_df = regressionResult(dataType; L0 = L0, f = f)

    p1 = plotEC50(L0, f; dataType = dataType, fit = res, Cellidx = 6, Recepidx = 2, Rbound = true)
    #p1 = mixEC50()
    draw(
        SVG("figureS1.svg", 10inch, 8inch),
        plotGrid((1,1), [p1]; widths = [3]),
    )
end
