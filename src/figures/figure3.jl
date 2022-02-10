""" Figure 3: deconvolve the receptor and cell type functionality """

function oneCellTypeOnlyR2(dataType; L0 = 1e-9, f = 4, murine = true, cellTypes = ["ncMO", "cMO"])
    # Only testing on exponential method
    df = murine ? importDepletion(dataType) : importHumanized(dataType)
    Xdf = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = cellTypes)
    res = fitRegNNLS(Xdf; murine = murine, cellTypes = cellTypes, exp_method = true)
    loo_res = regLOO(Xdf; murine = murine, cellTypes = cellTypes, exp_method = true)
    return res.R2, [l.R2 for l in loo_res]
end

function oneCellTypeOnlyplot(dataType; L0 = 1e-9, f = 4, murine = true)
    allCellTypes = murine ? murineCellTypes : humanCellTypes
    mR2, lR2 = oneCellTypeOnlyR2(dataType; L0 = L0, f = f, murine = murine, cellTypes = nothing)
    mR2s, lR2s = [mR2], [lR2]
    for ct in allCellTypes
        mR2, lR2 = oneCellTypeOnlyR2(dataType; L0 = L0, f = f, murine = murine, cellTypes = [ct])
        append!(mR2s, [mR2])
        append!(lR2s, [lR2])
    end
    looR2 = hcat(lR2s...)
    low = [lower(x) for x in eachslice(looR2, dims = 2)]
    high = [upper(x) for x in eachslice(looR2, dims = 2)]
    return plot(
        DataFrame(CellTypes = vcat(["All"], allCellTypes .* " only"), R2 = mR2s, ymin = low, ymax = high),
        x = "CellTypes",
        y = "R2",
        ymin = "ymin",
        ymax = "ymax",
        Geom.errorbar,
        Geom.bar,
        Guide.title("Regression <i>R</i><sup>2</sup> with single cell type"),
        Guide.xticks(orientation = :vertical),
        Guide.xlabel("Cell types"),
        Guide.ylabel("<i>R</i><sup>2</sup>"),
        style(bar_spacing = 10pt, stroke_color = c -> "black"),
    )
end

function figure3()
    setGadflyTheme()

    mp1, mp2, mp3 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, legend = true)
    mp4 = oneCellTypeOnlyplot("melanoma"; L0 = 1e-9, f = 6, murine = true)

    ip1, ip2, ip3 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, legend = true)
    ip4 = oneCellTypeOnlyplot("ITP"; L0 = 1e-8, f = 10, murine = true)

    draw(SVG("figure3.svg", 1600px, 600px), plotGrid((2, 4), [mp1, mp2, mp3, mp4, ip1, ip2, ip3, ip4]))
end
