""" Figure 1: show the mixture IC binding data """

function figureS1()
    setGadflyTheme()
    draw(PDF("figureS1.pdf", 20inch, 16inch), plotMixSubplots(x -> splot_origData(x; match_y = false), averageMixData()))
end
