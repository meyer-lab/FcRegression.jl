""" Figure A4: deconvolve the receptor and cell type functionality """

function oneCellTypeOnlyR2(dataType; L0 = 1e-9, f = 4, murine = true, cellTypes = ["ncMO", "cMO"])
    # Only testing on exponential method
    df = murine ? importDepletion(dataType) : importHumanized(dataType)
    Xdf = modelPred(df; L0 = L0, f = f, murine = murine, cellTypes = cellTypes)
    res = fitRegNNLS(Xdf; murine = murine, cellTypes = cellTypes)
    loo_res = regLOO(Xdf; murine = murine, cellTypes = cellTypes)
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

function figureA4()
    setGadflyTheme()

    mp1, mp2, mp3 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, legend = true)
    mp4 = oneCellTypeOnlyplot("melanoma"; L0 = 1e-9, f = 6, murine = true)

    ip1, ip2, ip3 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, legend = true)
    ip4 = oneCellTypeOnlyplot("ITP"; L0 = 1e-8, f = 10, murine = true)

    draw(PDF("figure4.pdf", 1600px, 600px), plotGrid((2, 4), [mp1, mp2, mp3, mp4, ip1, ip2, ip3, ip4]))
end


function figure4Fig5c()
    df = FcRegression.importHumanized("bloodFig5c")

    Kav0 = FcRegression.importKav(; murine = false)
    Kav1 = FcRegression.extractNewHumanKav()

    # validate
    opt0, _, _ = FcRegression.runRegMAP(df; murine = false, Kav = Kav0)
    c0, _ = FcRegression.runRegMCMC(df; murine = false, Kav = Kav0, mcmc_iter = 200)

    plmap0 = FcRegression.plotRegMCMC(opt0, deepcopy(df); ptitle = "Fig5c, MAP, old affinity", Kav = Kav0)
    plmc0 = FcRegression.plotRegMCMC(c0, deepcopy(df); ptitle = "Fig5c, MCMC, old affinity", Kav = Kav0)
    cpl0 = FcRegression.plotRegParams(c0; ptitle = "Fig5c, MCMC, old affinity", legend = true, Kav = Kav0)


    # new
    opt1, _, _ = FcRegression.runRegMAP(df; murine = false, Kav = Kav1)
    c1, _ = FcRegression.runRegMCMC(df; murine = false, Kav = Kav1, mcmc_iter = 200)

    plmap1 = FcRegression.plotRegMCMC(opt1, deepcopy(df); ptitle = "Fig5c, MAP, new affinity", Kav = Kav1)
    plmc1 = FcRegression.plotRegMCMC(c1, deepcopy(df); ptitle = "Fig5c, MCMC, new affinity", Kav = Kav1)
    cpl1 = FcRegression.plotRegParams(c1; ptitle = "Fig5c, MCMC, new affinity", legend = true, Kav = Kav1)



end
