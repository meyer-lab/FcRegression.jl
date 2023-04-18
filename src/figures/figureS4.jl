""" Supplemental figure to Fig. 5 """


function figureS4()
    setGadflyTheme()

    p0 = FcRegression.plotEffectorPredict(
        FcRegression.importEffectorBind(; avg = true),
        FcRegression.predictLbound(FcRegression.importKav(; murine = false)); 
        title = "Leukocyte binding with documented affinities"
    )

    p1 = FcRegression.plotEffectorPredict(
        FcRegression.importEffectorBind(; avg = true),
        FcRegression.predictLbound(FcRegression.extractNewHumanKav()); 
        title = "Leukocyte binding with updated affinities"
    )

    # FcgRIIIB: IgG2 ~= 60000, IgG4 ~= 80000

    c = FcRegression.inferCD16b(; mcmc_iter = 5000)
    af1, af2 = FcRegression.plotCD16bCHOaff(c[4000:end])

    draw(
        PDF("output/figureS4.pdf", 14inch, 3inch),
        plotGrid(
            (1, 4),
            [p0, p1, af1, af2];
        )
    )
end