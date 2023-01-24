""" Plot in vivo regression results. Return predictions, cell weights, and receptor weights for MAP and MCMC. """
function plot_regressions(df; Kav::DataFrame, murine = false, cellTypes = nothing, ptitle = "")
    opt, optcv, cdf = runRegMAP(df; murine = murine, Kav = Kav, cellTypes = cellTypes)
    c0, ccdf0 = FcRegression.runRegMCMC(df; murine = murine, Kav = Kav0, mcmc_iter = 200, cellTypes = cellTypes)

    pl_map = plotRegMCMC(opt, deepcopy(df); ptitle = ptitle * "[MAP]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    pl_mcmc = plotRegMCMC(c, deepcopy(df); ptitle = ptitle * "[MCMC]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    cell_map, act_map = plotRegParams(optcv; ptitle = ptitle * "[MAP]", legend = true, Kav = Kav, cellTypes = cellTypes)
    cell_mcmc, act_mcmc = plotRegParams(c; ptitle = ptitle * "[MCMC]", legend = true, Kav = Kav, cellTypes = cellTypes)
    return [pl_map, cell_map, act_map, pl_mcmc, cell_mcmc, act_mcmc]
end


function figure6(ssize = (8.5inch, 5.5inch); cellTypes = ["ncMO", "cMO", "Neu"], mcmc_iter = 50000, suffix = "0117_", kwargs...)
    setGadflyTheme()
    df = FcRegression.importHumanized("ITP")

    Kav0 = FcRegression.extractNewHumanKav(; old = true)
    Kav1 = FcRegression.extractNewHumanKav(; old = false)

    c1, ccdf1 = FcRegression.runRegMCMC(
        df,
        "regMCMC_$(suffix)1.dat";
        murine = false,
        Kav = Kav1,
        fitActI = false,
        mcmc_iter = mcmc_iter,
        cellTypes = cellTypes,
    )
    c0, ccdf0 = FcRegression.runRegMCMC(
        df,
        "regMCMC_$(suffix)0.dat";
        murine = false,
        Kav = Kav0,
        fitActI = false,
        mcmc_iter = mcmc_iter,
        cellTypes = cellTypes,
    )

    c0 = c0[(mcmc_iter ÷ 10 * 7):end]
    c1 = c1[(mcmc_iter ÷ 10 * 7):end]

    pl_map0 = FcRegression.plotRegMCMC(
        c0,
        deepcopy(df);
        ptitle = "documented affinities",
        colorL = "Genotype",
        shapeL = "Condition",
        legend = false,
        Kav = Kav0,
        cellTypes = cellTypes,
    )
    cell_map0, act_map0 =
        FcRegression.plotRegParams(c0; ptitle = "documented affinities", legend = false, Kav = Kav0, cellTypes = cellTypes, cell_max_y = 1.5)

    pl_map1 = FcRegression.plotRegMCMC(
        c1,
        deepcopy(df);
        ptitle = "updated affinities",
        colorL = "Genotype",
        shapeL = "Condition",
        legend = true,
        Kav = Kav1,
        cellTypes = cellTypes,
    )
    cell_map1, act_map1 =
        FcRegression.plotRegParams(c1; ptitle = "updated affinities", legend = true, Kav = Kav1, cellTypes = cellTypes, cell_max_y = 1.5)

    pl_mapL = FcRegression.plotRegMCMC(
        c1,
        deepcopy(df);
        ptitle = "updated affinities",
        colorL = "Genotype",
        shapeL = "Condition",
        legend = true,
        Kav = Kav1,
        cellTypes = cellTypes,
    )
    cell_mapL, act_mapL = FcRegression.plotRegParams(c1; ptitle = "updated affinities", legend = true, Kav = Kav1, cellTypes = cellTypes)

    pl = FcRegression.plotGrid(
        (2, 3),
        [nothing, pl_map0, pl_map1, nothing, cell_map0, cell_map1];
        sublabels = "abc de",
        widths = [1 1 1.35; 1 1 1.35],
        kwargs...,
    )
    draw(PDF("output/figure6_$suffix.pdf", ssize[1], ssize[2]), pl)
end
