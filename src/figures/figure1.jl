""" Figure 1: show the mixture IC binding data """

function figure1()
    setGadflyTheme()
    draw(SVG("figure1.svg", 2500px, 1000px), plotMixOriginalData())
end
