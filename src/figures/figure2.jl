function figure2()
    setGadflyTheme()
    fig_melanoma = figureW("melanoma"; IgGX = 2, IgGY = 3, L0 = 1e-9, f = 6, murine = true, legend = false)
    fig_ITP = figureW("ITP"; IgGX = 1, IgGY = 2, L0 = 1e-9, f = 6, murine = true, legend = true)

    draw(
        SVG("figure2.svg", 7inch, 7inch),
        plotGrid(
            (2, 2),
            [fig_melanoma[1] fig_ITP[1] fig_melanoma[3] fig_ITP[3]];
            widths = [3 4; 1 1],
        ),
    )

    draw(
        SVG("figureS2.svg", 9inch, 6inch),
        plotGrid(
            (2, 2),
            [fig_melanoma[2] fig_melanoma[5]; fig_ITP[2] fig_ITP[5]],
            widths = [3 4; 1 1],
            heights = [2, 1],
        ),
    )
end
