""" Figure 3: predicting in vivo effect """

function figure3()
    setGadflyTheme()

    df = MixtureCellSeparateFit(loadMixData(); logscale = false)
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    avp_mel, _, cell_mel, recep_mel = figureW("melanoma"; L0 = 1e-9, f = 6)
    avp_itp, _, cell_itp, recep_itp = figureW("ITP"; L0 = 1e-8, f = 10)

    draw(SVG("figure3.svg", 1300px, 600px), plotGrid((2, 4), [nothing, avp_mel, cell_mel, recep_mel, nothing, avp_itp, cell_itp, recep_itp]))
end