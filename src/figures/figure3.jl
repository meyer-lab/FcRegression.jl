""" Figure 3: deconvolve the receptor and cell type functionality """

function compareRegMethodPlot(dataType; L0, f, murine = true)
    res1, loo_res1, odf1 = regressionResult(dataType; L0 = L0, f = f, murine = murine, exp_method = false, fit_ActI = false)
    res2, loo_res2, odf2 = regressionResult(dataType; L0 = L0, f = f, murine = murine, exp_method = false, fit_ActI = true)
    res3, loo_res3, odf3 = regressionResult(dataType; L0 = L0, f = f, murine = murine, exp_method = true, fit_ActI = false)
    res4, loo_res4, odf4 = regressionResult(dataType; L0 = L0, f = f, murine = murine, exp_method = true, fit_ActI = true)

    # plot R2 for old and new methods
    loolist = [loo_res1, loo_res2, loo_res3, loo_res4]
    R2 = [lool.R2 for lool in [res1, res2, res3, res4]]
    R2_lower = [lower([loo.R2 for loo in lool]) for lool in loolist]
    R2_upper = [upper([loo.R2 for loo in lool]) for lool in loolist]
    R2df = DataFrame(Method = ["Original", "ActI fit", "Exp only", "Exp and ActI fit"], R2 = R2, ymin = R2_lower, ymax = R2_upper)
    plR2 = plot(R2df, 
        x = "Method", y = "R2", ymin = "ymin", ymax = "ymax", 
        Geom.errorbar,
        Stat.dodge(axis=:x),
        Geom.bar(position=:dodge),
        Scale.x_discrete(),
        Scale.y_continuous(minvalue = 0.0),
        Guide.title("Fitting R<sup>2</sup> for different methods, $dataType"),
        style(bar_spacing = 5mm, stroke_color=c->"black"),
    )

    supppls = []
    for rres in [[res1, loo_res1, odf1], [res2, loo_res2, odf2], [res3, loo_res3, odf3]]
        p1, p2, p3 = figureW(rres..., dataType; L0 = L0, f = f, murine = murine, legend = true)
        append!(supppls, [p1, p2, p3])
    end
    p1, p2, p3 = figureW(res4, loo_res4, odf4, dataType; L0 = L0, f = f, murine = murine, legend = true)
    
    return plR2, [p1, p2, p3], supppls
end

function figure3()
    setGadflyTheme()
    
    MplR2, Mpls, Msupp = compareRegMethodPlot("melanoma"; L0 = 1e-9, f = 6, murine = true)
    IplR2, Ipls, Isupp = compareRegMethodPlot("ITP"; L0 = 1e-8, f = 10, murine = true)
    
    draw(SVG("figure3.svg", 1300px, 600px), plotGrid((2, 4), [MplR2, Mpls..., IplR2, Ipls...]))
    draw(SVG("figureS3.svg", 900px, 1800px), plotGrid((6, 3), [Msupp..., Isupp...]))
end

function figureS3()
    setGadflyTheme()

    A1, C1, R1 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = false, fit_ActI = false)
    A2, C2, R2 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = false, fit_ActI = true)
    A3, C3, R3 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = true, fit_ActI = false)
    A4, C4, R4 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = true, fit_ActI = true)
    
    A01, C01, R01 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = false, fit_ActI = false)
    A02, C02, R02 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = false, fit_ActI = true)
    A03, C03, R03 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = true, fit_ActI = false)
    A04, C04, R04 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = true, fit_ActI = true)

    draw(SVG("figureS3.svg", 900px, 2400px), plotGrid((8, 3), [A1,C1,R1,A2,C2,R2,A3,C3,R3,A4,C4,R4,A01,C01,R01,A02,C02,R02,A03,C03,R03,A04,C04,R04]))
end