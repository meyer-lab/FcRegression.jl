function figureDS(dataType; L0 = 1e-9, f = 4, murine::Bool)

    if murine
        df = importDepletion(dataType)
    else
        df = importHumanized(dataType)
    end
    
    res, odf, effects, ActI_df = regressionResult(dataType; L0 = L0, f = f, murine = murine)

    setGadflyTheme()
    p1 = plotDepletionSynergy(
        1,
        2;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p2 = plotDepletionSynergy(
        1,
        3;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p3 = plotDepletionSynergy(
        1,
        4;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p4 = plotDepletionSynergy(
        1,
        5;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p5 = plotDepletionSynergy(
        2,
        3;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p6 = plotDepletionSynergy(
        2,
        4;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p7 = plotDepletionSynergy(
        2,
        5;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p8 = plotDepletionSynergy(
        3,
        4;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p9 = plotDepletionSynergy(
        3,
        5;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )
    p10 = plotDepletionSynergy(
        4,
        5;
        L0 = L0,
        f = f,
        murine = murine,
        fit = res,
        neutralization = ("Neutralization" in names(df)),
        c1q = ("C1q" in effects.Component),
        dataType = dataType,
    )

    return p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
end

function figure_MSyn_melanoma()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureDS("melanoma"; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_Mmelanoma.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_MSyn_ITP()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureDS("ITP"; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_MITP.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end
