function figure6()
    df = FcRegression.importHumanized("ITP")

    Kav0 = FcRegression.importKav(; murine = false)
    Kav1 = FcRegression.extractNewHumanKav()

    #FcRegression.modelPred(df, murine = false, Kav = Kav0)

    opt0, _, _ = FcRegression.runRegMAP(df; murine = false, Kav = Kav0)

    plmap0 = FcRegression.plotRegMCMC(opt0, deepcopy(df); ptitle = "Schwab, MAP, old affinity", Kav = Kav0)


    # validate
    opt0, opt0l, _ = FcRegression.runRegMAP(df; murine = false, Kav = Kav0)
    c0, _ = FcRegression.runRegMCMC(df; murine = false, Kav = Kav0, mcmc_iter = 200)

    plmap0 = FcRegression.plotRegMCMC(opt0, deepcopy(df); ptitle = "Schwab, MAP, old affinity", 
        Kav = Kav0, colorL="Genotype", shapeL="Condition")
    plmc0 = FcRegression.plotRegMCMC(c0, deepcopy(df); ptitle = "Schwab, MCMC, old affinity", 
        Kav = Kav0, colorL="Genotype", shapeL="Condition")
    cpl0 = FcRegression.plotRegParams(c0; ptitle = "Schwab, MCMC, old affinity", legend = true, Kav = Kav0)


    # new
    opt1, opt1l, _ = FcRegression.runRegMAP(df; murine = false, Kav = Kav1)
    c1, _ = FcRegression.runRegMCMC(df; murine = false, Kav = Kav1, mcmc_iter = 200)

    plmap1 = FcRegression.plotRegMCMC(opt1, deepcopy(df); ptitle = "Schwab, MAP, new affinity",
        Kav = Kav1, colorL="Genotype", shapeL="Condition")
    plmc1 = FcRegression.plotRegMCMC(c1, deepcopy(df); ptitle = "Schwab, MCMC, new affinity", Kav = Kav1)
    cpl1 = FcRegression.plotRegParams(c1; ptitle = "Schwab, MCMC, new affinity", legend = true, Kav = Kav1)


end