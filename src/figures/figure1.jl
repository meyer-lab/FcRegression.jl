""" Figure 1: show the mixture IC binding data """

function figure1()
    setGadflyTheme()

    df = MixtureCellSeparateFit(loadMixData(); logscale = false)
    draw(SVG("figure1.svg", 2500px, 1000px), makeMixturePairSubPlots(df; logscale = false))
end
