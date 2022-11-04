""" Figure S2: Compare tanh and exp regression """

function figureS3()
    df = FcRegression.importDepletion("melanoma")

    Kav = FcRegression.importKav(; murine = true, retdf = true, IgG2bFucose = true)
    Kav = Kav[.!(contains.(Kav."IgG", "SA")), :]

    c0, ccdf0 = FcRegression.runRegMCMC(df; 
        murine = true, Kav = Kav, fitActI = false, mcmc_iter = 100_000, link = FcRegression.exponential, cellTypes = ["ncMO", "cMO", "NKs", "Neu", "EO"])

    c1, ccdf1 = FcRegression.runRegMCMC(df; 
        murine = true, Kav = Kav, fitActI = false, mcmc_iter = 100_000, link = FcRegression.tanh, cellTypes = ["ncMO", "cMO", "NKs", "Neu", "EO"])
    
    pl_map0 = FcRegression.plotRegMCMC(
        c0,
        deepcopy(df);
        ptitle = "Exponential",
        colorL = "Background",
        shapeL = "Condition",
        legend = true,
        Kav = Kav,
        link = FcRegression.exponential,
    )

    pl_map1 = FcRegression.plotRegMCMC(
        c1,
        deepcopy(df);
        ptitle = "tanh",
        colorL = "Background",
        shapeL = "Condition",
        legend = true,
        Kav = Kav,
        link = FcRegression.tanh,
    )

    pl = FcRegression.plotGrid((1, 2), [pl_map1, pl_map0];)
    draw(PDF("output/figureS3.pdf", 7inch, 3inch), pl)
end