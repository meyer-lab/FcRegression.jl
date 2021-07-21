""" Figure 3: deconvolve the receptor and cell type functionality """

function figure3()
    setGadflyTheme()

    df = MixtureCellSeparateFit(loadMixData(); logscale = false)
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    avp_mel, _, cell_mel, recep_mel, _, _ = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true)
    avp_itp, _, cell_itp, recep_itp, _, _ = figureW("ITP"; L0 = 1e-8, f = 10, murine = true)

    draw(SVG("figure3.svg", 1300px, 600px), plotGrid((2, 4), [nothing, avp_mel, cell_mel, recep_mel, nothing, avp_itp, cell_itp, recep_itp]))
end
