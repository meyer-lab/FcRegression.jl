""" Figure 1: show the mixture IC binding data """

function figureS1()
    setGadflyTheme()
    draw(SVG("figureS1.svg", 20inch, 16inch), plotMixSubplots(splot_origData, loadMixData(); avg = true))
    draw(SVG("figureS1_pred.svg", 25inch, 16inch), plotMixSubplots(splot_contPred, loadMixData(); logscale = false))
end
