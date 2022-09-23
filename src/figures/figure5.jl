""" Plot in vivo regression results. Return predictions, cell weights, and receptor weights for MAP and MCMC. """
function plot_regressions(df; Kav::DataFrame, murine = false, cellTypes = nothing, ptitle = "")
    opt, optcv, cdf = runRegMAP(df; murine = murine, Kav = Kav, cellTypes = cellTypes)
    c, ccdf = runRegMCMC(df; murine = murine, Kav = Kav, mcmc_iter = 200, cellTypes = cellTypes)

    pl_map = plotRegMCMC(opt, deepcopy(df); ptitle = ptitle * "[MAP]", Kav = Kav, 
        cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    pl_mcmc = plotRegMCMC(c, deepcopy(df); ptitle = ptitle * "[MCMC]", Kav = Kav, cellTypes = cellTypes, colorL = "Genotype", shapeL = "Condition")
    cell_map, act_map = plotRegParams(optcv; ptitle = ptitle * "[MAP]", legend = true, Kav = Kav, cellTypes = cellTypes)
    cell_mcmc, act_mcmc = plotRegParams(c; ptitle = ptitle * "[MCMC]", legend = true, Kav = Kav, cellTypes = cellTypes)
    return [pl_map, cell_map, act_map, pl_mcmc, cell_mcmc, act_mcmc]
end


function figure5()
    df = FcRegression.importHumanized("ITP");

    Kav0 = FcRegression.importKav(; murine = false);
    Kav1 = FcRegression.extractNewHumanKav();
    Kav0 = Kav0[!, Not(["FcgRIIB-232T", "FcgRIIC-13N"])];
    Kav1 = Kav1[!, Not(["FcgRIIB-232T", "FcgRIIC-13N"])];

    opt0, optcv0, mapdf0 = FcRegression.runRegMAP(df, "MAP_0923_A0.dat"; murine = false, Kav = Kav0);
    opt1, optcv1, mapdf1 = FcRegression.runRegMAP(df, "MAP_0923_A1.dat"; murine = false, Kav = Kav1);

    pl_map0 = FcRegression.plotRegMCMC(opt0, deepcopy(df); ptitle = "documented affinities", 
        Kav = Kav0, colorL = "Genotype", shapeL = "Condition", legend = false);
    cell_map0, act_map0 = FcRegression.plotRegParams(optcv0; ptitle = "documented affinities", legend = false, Kav = Kav0);

    pl_map1 = FcRegression.plotRegMCMC(opt1, deepcopy(df); ptitle = "updated affinities", 
        Kav = Kav1, colorL = "Genotype", shapeL = "Condition", legend = false);
    cell_map1, act_map1 = FcRegression.plotRegParams(optcv1; ptitle = "updated affinities", legend = false, Kav = Kav1);


    pl = FcRegression.plotGrid((2, 4), [nothing, pl_map0, cell_map0, act_map0, nothing, pl_map1, cell_map1, act_map1]; 
        sublabels = "abdf ceg")
    return draw(PDF("figure5.pdf", 12inch, 6inch), pl)
end
