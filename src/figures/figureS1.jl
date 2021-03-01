""" Figure S1 """

function figureS1()

    setGadflyTheme()

    dataType = "ITP"
    df = importDepletion(dataType)
    L0 = 1e-9
    f = 4
    murine = true

    res, odf, Cell_df, ActI_df = regressionResult(dataType; L0 = L0, f = f, murine = murine)

    p1 = plotEC50(L0, f, 6, 2; dataType = dataType, fit = res, Rbound = true)
    p2 = maxSynergyBar(L0, f; dataType = dataType, fit = res, Rbound = false)
    #p1 = mixEC50()
    draw(SVG("figureS1.svg", 5inch, 4inch), plotGrid((1, 2), [p1 p2]; widths = [3 3]))
end
