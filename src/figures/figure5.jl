""" Figure 5: validating the in vivo effect of mixtures """

function figure5()
    setGadflyTheme()
    draw(
        SVG("figure5.svg", 2inch, 2inch),
        context(),
    )
end
