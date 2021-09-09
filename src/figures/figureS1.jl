""" Figure 1: show the mixture IC binding data """

function figureS1()
    setGadflyTheme()
    draw(SVG("figureS1.svg", 20inch, 16inch), plotMixSubplots(splot_origData, averageMixData()))
    draw(SVG("figureS3.svg", 20inch, 16inch), plotMixSubplots(splot_contPred, loadMixData()))
end
