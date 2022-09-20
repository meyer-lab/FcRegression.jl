""" Plot in vivo regression results. Return predictions, cell weights, and receptor weights for MAP and MCMC. """
function plot_regressions(df; Kav::DataFrame, murine = false, cellTypes = nothing, ptitle = "")
    opt, optcv, cdf = runRegMAP(df; murine = murine, Kav = Kav, cellTypes = cellTypes)
    c, ccdf = runRegMCMC(df; murine = murine, Kav = Kav, mcmc_iter = 200, cellTypes = cellTypes)

    pl_map = plotRegMCMC(opt, deepcopy(df); ptitle = ptitle * "[MAP]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    pl_mcmc = plotRegMCMC(c, deepcopy(df); ptitle = ptitle * "[MCMC]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    cell_map, act_map = plotRegParams(optcv; ptitle = ptitle * "[MAP]", legend = true, Kav = Kav, cellTypes = cellTypes)
    cell_mcmc, act_mcmc = plotRegParams(c; ptitle = ptitle * "[MCMC]", legend = true, Kav = Kav, cellTypes = cellTypes)
    return [pl_map, cell_map, act_map, pl_mcmc, cell_mcmc, act_mcmc]
end


function figure6()
    df = FcRegression.importHumanized("ITP")

    Kav0 = FcRegression.importKav(; murine = false)
    Kav1 = FcRegression.extractNewHumanKav()

    pls0a = FcRegression.plot_regressions(df; Kav = Kav1, ptitle = "Schwab, old Kav(A)")
    pls1a = FcRegression.plot_regressions(df; Kav = Kav1, ptitle = "Schwab, new Kav(A)")
    pls0t = FcRegression.plot_regressions(df; Kav = Kav0, ptitle = "Schwab, old Kav(3)")
    pls1t = FcRegression.plot_regressions(df; Kav = Kav1, ptitle = "Schwab, new Kav(3)")

    draw(PDF("old Kav, all cells.pdf", 12inch, 8inch), FcRegression.plotGrid((2, 3), pls0a))
    draw(PDF("new Kav, all cells.pdf", 12inch, 8inch), FcRegression.plotGrid((2, 3), pls1a))
    draw(PDF("old Kav, three cells.pdf", 12inch, 8inch), FcRegression.plotGrid((2, 3), pls0t))
    draw(PDF("new Kav, three cells.pdf", 12inch, 8inch), FcRegression.plotGrid((2, 3), pls1t))
end
