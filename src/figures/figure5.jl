""" Figure 5: validating the in vivo effect of mixtures """

function plot_mix_ITP_depletion()
    df = invivo_prediction(importDeplExp(); L0 = 1e-9, f = 4, KxStar = KxConst)
    df[!, "Dose"] = [generate_mixture_name(dfr) for dfr in eachrow(df)]
    R2anno = "<i>R</i><sup>2</sup>" * @sprintf("=%.3f", R2(df."depletion", df."Predicted"; logscale = false))
    setGadflyTheme()
    return plot(
        df,
        x = "depletion",
        y = "Predicted",
        color = "Dose",
        Geom.point,
        Geom.abline(color = "green"),
        Guide.title("Mice <i>in vivo</i> platelet depletion"),
        Guide.xlabel("Measured"),
        Guide.annotation(compose(context(), text(-0.6, 0.7, R2anno), fontsize(12pt), font("Helvetica"))),
    )
end

function figure5()
    setGadflyTheme()
    pl = plot_mix_ITP_depletion()
    draw(SVG("figure5.svg", 6inch, 4inch), pl)
end
