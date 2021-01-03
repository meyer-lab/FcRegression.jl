function figure4()

    setGadflyTheme()
    ITP_Dep = figureW("ITP"; IgGX = 2, IgGY = 4, L0 = 1e-9, murine = true, legend = true)
    ITP_Kupffer_Act = figureW("ITP"; IgGX = 2, IgGY = 4, L0 = 1e-9, murine = true, legend = true, Cellidx = 6)
    ITP_Kupffer_Bound = figureW("ITP"; IgGX = 2, IgGY = 4, L0 = 1e-9, murine = true, legend = true, Cellidx = 6, Rbound = true)

    draw(
        SVG("figure4.svg", 10inch, 8inch),
        plotGrid((2, 3), [ITP_Dep[4] ITP_Kupffer_Act[4] ITP_Kupffer_Bound[4] ITP_Dep[6] ITP_Kupffer_Act[6] ITP_Kupffer_Bound[6]]; widths = [3 3 3; 4 4 4]),
    )
end
