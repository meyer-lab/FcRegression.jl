""" Figure 3: deconvolve the receptor and cell type functionality """

function figure3()
    setGadflyTheme()

    avp_mel, _, cell_mel, recep_mel = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true)
    avp_itp, _, cell_itp, recep_itp = figureW("ITP"; L0 = 1e-8, f = 10, murine = true)

    draw(SVG("figure3.svg", 1300px, 600px), plotGrid((2, 4), [nothing, avp_mel, cell_mel, recep_mel, nothing, avp_itp, cell_itp, recep_itp]))
end

function figureS3()
    setGadflyTheme()

    A1, _, C1, R1 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = false, fit_ActI = false)
    A2, _, C2, R2 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = false, fit_ActI = true)
    A3, _, C3, R3 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = true, fit_ActI = false)
    A4, _, C4, R4 = figureW("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = true, fit_ActI = true)
    
    A01, _, C01, R01 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = false, fit_ActI = false)
    A02, _, C02, R02 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = false, fit_ActI = true)
    A03, _, C03, R03 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = true, fit_ActI = false)
    A04, _, C04, R04 = figureW("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = true, fit_ActI = true)

    draw(SVG("figureS3.svg", 900px, 2400px), plotGrid((8, 3), [A1,C1,R1,A2,C2,R2,A3,C3,R3,A4,C4,R4,A01,C01,R01,A02,C02,R02,A03,C03,R03,A04,C04,R04]))
end