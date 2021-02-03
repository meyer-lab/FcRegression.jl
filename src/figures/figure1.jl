""" Figure 1: show the mixture IC binding data """

function figure1()
    setGadflyTheme()
    draw(SVG("figure1.svg", 2500px, 1000px), plotMixOriginalData(loadMixData()))
    draw(SVG("figure1A.svg", 2500px, 1000px), plotMixOriginalData(PCAData(; cutoff = 0.95)))
    draw(SVG("figure1B.svg", 2500px, 1000px), plotMixOriginalData(PCAData(; cutoff = 0.9)))
    draw(SVG("figure1C.svg", 2500px, 1000px), plotMixOriginalData(PCAData(; cutoff = 0.8)))
end


""
