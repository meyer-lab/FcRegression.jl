""" Supplemental figure to Fig. 5 """


function figureS4()

    df = FcRegression.importEffectorBind(; avg = true)

    # TODO: why no errorbar?

    
    Kav0 = FcRegression.importKav(; murine = false);

    Kav1 = FcRegression.extractNewHumanKav();

    plotEffectorPredict(df, predictLbound(Kav0); title = nothing)
    
    plotEffectorPredict(df, predictLbound(Kav1); title = nothing)


    FcRegression.plotEffectorPredict(
        FcRegression.importEffectorBind(; avg = true),
        FcRegression.predictLbound(); 
        title = nothing
    )

    # FcgRIIIB: IgG2 ~= 60000, IgG4 ~= 80000



end