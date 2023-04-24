""" Plot in vivo regression results. Return predictions, cell weights, and receptor weights for MAP and MCMC. """
function plot_regressions(df; Kav::DataFrame, murine = false, cellTypes = nothing, ptitle = "")
    opt, optcv, cdf = runRegMAP(df; murine = murine, Kav = Kav, cellTypes = cellTypes)
    c0, ccdf0 = runRegMCMC(df; murine = murine, Kav = Kav0, mcmc_iter = 200, cellTypes = cellTypes)

    pl_map = plotRegMCMC(opt, deepcopy(df); ptitle = ptitle * "[MAP]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    pl_mcmc = plotRegMCMC(c, deepcopy(df); ptitle = ptitle * "[MCMC]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    cell_map, act_map = plotRegParams(optcv; ptitle = ptitle * "[MAP]", legend = true, Kav = Kav, cellTypes = cellTypes)
    cell_mcmc, act_mcmc = plotRegParams(c; ptitle = ptitle * "[MCMC]", legend = true, Kav = Kav, cellTypes = cellTypes)
    return [pl_map, cell_map, act_map, pl_mcmc, cell_mcmc, act_mcmc]
end


function figure6(ssize = (8.5inch, 5.5inch); cellTypes = ["ncMO", "cMO", "Neu"], mcmc_iter = 5000, suffix = "0117_", kwargs...)
    setGadflyTheme()
    df = importHumanized()

    Kav0 = importKavDist(; murine = false, regularKav = true, retdf = true)
    Kav1 = extractNewHumanKav()

    c1, ccdf1 = runRegMCMC(
        df,
        "regMCMC_$(suffix)1.dat";
        murine = false,
        Kav = Kav1,
        fitActI = false,
        mcmc_iter = mcmc_iter,
        cellTypes = cellTypes,
    )
    c0, ccdf0 = runRegMCMC(
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

    pl_map0 = plotRegMCMC(
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
        plotRegParams(c0; ptitle = "documented affinities", legend = false, Kav = Kav0, cellTypes = cellTypes, cell_max_y = 1.5)

    pl_map1 = plotRegMCMC(
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
        plotRegParams(c1; ptitle = "updated affinities", legend = true, Kav = Kav1, cellTypes = cellTypes, cell_max_y = 1.5)

    pl_mapL = plotRegMCMC(
        c1,
        deepcopy(df);
        ptitle = "updated affinities",
        colorL = "Genotype",
        shapeL = "Condition",
        legend = true,
        Kav = Kav1,
        cellTypes = cellTypes,
    )
    cell_mapL, act_mapL = plotRegParams(c1; ptitle = "updated affinities", legend = true, Kav = Kav1, cellTypes = cellTypes)

    pl = plotGrid(
        (2, 3),
        [nothing, pl_map0, pl_map1, nothing, cell_map0, cell_map1];
        sublabels = "abc de",
        widths = [1 1 1.35; 1 1 1.35],
        kwargs...,
    )
    draw(PDF("output/figure6_$suffix.pdf", ssize[1], ssize[2]), pl)
end
