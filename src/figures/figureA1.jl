using ColorSchemes




function figureA1()
    setGadflyTheme()

    p1, p2 = bindVSaff()

    # Specific IgG pair - receptor interaction
    df = averageMixData()
    igg12_1 = splot_origData(df[(df."Receptor" .== "FcgRI") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG2"), :]; match_y = false)
    igg14_1 =
        splot_origData(df[(df."Receptor" .== "FcgRIIIA-158F") .& (df."subclass_1" .== "IgG1") .& (df."subclass_2" .== "IgG4"), :]; match_y = false)

    score, loading, vars_expl = mixtureDataPCA()
    vars = plot(
        DataFrame(Components = 1:length(vars_expl), R2X = vars_expl),
        x = "Components",
        y = "R2X",
        label = [@sprintf("%.2f %%", i * 100) for i in vars_expl],
        Geom.point,
        Geom.line,
        Geom.label,
        Scale.x_discrete,
        Scale.y_continuous(minvalue = 0.5),
        Guide.title("Variance Explained by PCA"),
        Guide.xlabel("Number of components"),
        Guide.ylabel("R2X"),
    )

    SP4 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs", xx = "PC 1", yy = "PC 2")
    SP33 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs", xx = "PC 1", yy = "PC 2")
    SP4_13 = plot_PCA_score(score[score."Valency" .== 4, :]; title = "PCA Score, 4-valent ICs", xx = "PC 1", yy = "PC 3")
    SP33_13 = plot_PCA_score(score[score."Valency" .== 33, :]; title = "PCA Score, 33-valent ICs", xx = "PC 1", yy = "PC 3")
    LP = plot(
        loading,
        x = "PC 1",
        y = "PC 2",
        color = "Receptor",
        label = "Receptor",
        Geom.point,
        Geom.label,
        Guide.title("PCA Loadings"),
        Scale.x_continuous(minvalue = -1.0, maxvalue = 1.0),
        Scale.y_continuous(minvalue = -1.0, maxvalue = 1.0),
    )
    LP_13 = plot(
        loading,
        x = "PC 1",
        y = "PC 3",
        color = "Receptor",
        label = "Receptor",
        Geom.point,
        Geom.label,
        Guide.title("PCA Loadings"),
        Scale.x_continuous(minvalue = -1.0, maxvalue = 1.0),
        Scale.y_continuous(minvalue = -1.0, maxvalue = 1.0),
    )

    pl = plotGrid((3, 5), [nothing, nothing, p1, p2, vars, SP4, SP33, LP, igg12_1, igg14_1, SP4_13, SP33_13, LP_13];)
    return draw(PDF("figure1.pdf", 20inch, 10inch), pl)
end
