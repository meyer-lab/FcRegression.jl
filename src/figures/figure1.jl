""" Figure 1: show the mixture IC binding data """

function figure1()
    setGadflyTheme()
    draw(SVG("figure1.svg", 20inch, 16inch), plotMixOriginalData(loadMixData()))
    #draw(SVG("figure1A.svg", 20inch, 16inch), plotMixOriginalData(PCAData(; cutoff = 0.95)))
end


""
