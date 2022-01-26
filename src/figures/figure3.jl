""" Figure 3: deconvolve the receptor and cell type functionality """

function oneCellTypeOnlyR2(dataType; L0 = 1e-9, f = 4, murine = true, cellTypes = ["ncMO", "cMO"])
    df = murine ? importDepletion(dataType) : importHumanized(dataType)
    Xfc, Xdf, Y = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = cellTypes)
    res = fitRegression(Xfc, Xdf, Y; murine = murine, exp_method = true)
    loo_res = LOOCrossVal(Xfc, Xdf, Y; murine = murine)
    rY = exponential(regressionPred(Xfc, Xdf, res; murine = murine, cellTypes = cellTypes))
    loo_rYs = [exponential(regressionPred(Xfc, Xdf, rres; murine = murine, cellTypes = cellTypes)) for rres in loo_res]
    loo_R2s = [R2(Y, loo_rY; logscale = false) for loo_rY in loo_rYs]
    return R2(Y, rY; logscale = false) #, loo_R2s
end

function oneCellTypeOnlyplot(dataType; L0 = 1e-9, f = 4, murine = true)
    allCellTypes = murine ? murineCellTypes : humanCellTypes
    R2s = [oneCellTypeOnlyR2(dataType; L0 = L0, f = f, murine = murine, cellTypes = nothing)]
    for ct in allCellTypes
        append!(R2s, [oneCellTypeOnlyR2(dataType; L0 = L0, f = f, murine = murine, cellTypes = [ct])])
    end
    return plot(
            DataFrame(CellTypes=vcat(["All"], allCellTypes .* " only"), R2=R2s),
            x = "CellTypes",
            y = "R2",
            Geom.bar,
            Guide.title("Regression R<sup>2</sup> with single cell type"),
            Guide.xticks(orientation=:vertical),
            Guide.xlabel("Cell types"),
            Guide.ylabel("R<sup>2</sup>"),
            style(bar_spacing = 5mm),
        )
end

function figure3()
    setGadflyTheme()

    mres, mloo_res, modf = regressionResult("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = true, fit_ActI = true)
    mp1, mp2, mp3 = figureW(mres, mloo_res, modf, "melanoma"; L0 = 1e-9, f = 6, murine = true, legend = true)
    mp4 = oneCellTypeOnlyplot("melanoma"; L0 = 1e-9, f = 6, murine = true)

    ires, iloo_res, iodf = regressionResult("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = true, fit_ActI = true)
    ip1, ip2, ip3 = figureW(ires, iloo_res, iodf, "ITP"; L0 = 1e-8, f = 10, murine = true, legend = true)
    ip4 = oneCellTypeOnlyplot("ITP"; L0 = 1e-8, f = 10, murine = true)

    draw(SVG("figure3.svg", 1600px, 600px), plotGrid((2, 5), [nothing, mp1, mp2, mp3, mp4, nothing, ip1, ip2, ip3, ip4]))
end
