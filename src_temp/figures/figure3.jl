function figure3()
    setGadflyTheme()
    fig_ITP = figureW("ITP"; IgGX = 1, IgGY = 2, L0 = 1e-9, f = 6, murine = true, legend = true)

    draw(SVG("figure3.svg", 7inch, 7inch), plotGrid((2, 2), [fig_ITP[4] fig_ITP[6]];))
end
