function figure1()
    setGadflyTheme()
    df3 = PairFit(loadMixData(); logscale = false)

    draw(SVG("figure1.svg", 2500px, 1000px), makeMixturePairSubPlots(df3; logscale = false))
end
