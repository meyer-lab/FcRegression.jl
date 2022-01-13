""" Figure 3: deconvolve the receptor and cell type functionality """

function figure3()
    setGadflyTheme()

    mres, mloo_res, modf = regressionResult("melanoma"; L0 = 1e-9, f = 6, murine = true, exp_method = true, fit_ActI = true)
    mp1, mp2, mp3 = figureW(mres, mloo_res, modf, "melanoma"; L0 = 1e-9, f = 6, murine = true, legend = true)

    ires, iloo_res, iodf = regressionResult("ITP"; L0 = 1e-8, f = 10, murine = true, exp_method = true, fit_ActI = true)
    ip1, ip2, ip3 = figureW(ires, iloo_res, iodf, "ITP"; L0 = 1e-8, f = 10, murine = true, legend = true)

    draw(SVG("figure3.svg", 1300px, 600px), plotGrid((2, 4), [nothing, mp1, mp2, mp3, nothing, ip1, ip2, ip3]))
end
