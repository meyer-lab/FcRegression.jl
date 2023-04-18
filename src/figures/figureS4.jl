""" Supplemental figure to Fig. 5 """


function figureS4()
    setGadflyTheme()

    p0 = FcRegression.plotEffectorPredict(
        FcRegression.importEffectorBind(; avg = true),
        FcRegression.predictLbound(FcRegression.importKav(; murine = false)); 
        title = nothing
    )

    p1 = FcRegression.plotEffectorPredict(
        FcRegression.importEffectorBind(; avg = true),
        FcRegression.predictLbound(FcRegression.extractNewHumanKav()); 
        title = nothing
    )

    # FcgRIIIB: IgG2 ~= 60000, IgG4 ~= 80000

    draw(
        PDF("output/figureS4.pdf", 7inch, 3inch),
        plotGrid(
            (1, 2),
            [p0, p1];
        )
    )
end