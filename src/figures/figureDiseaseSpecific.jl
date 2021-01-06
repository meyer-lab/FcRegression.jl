function figure_melanoma()
    p1, p2, p3, p4, p5, p6 = figureW("melanoma"; IgGX = 3, IgGY = 5, L0 = 1e-9, f = 6, murine = true)
    draw(SVG("figure_Mmelanoma.svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_ITP()
    p1, p2, p3, p4, p5, p6 = figureW("ITP"; IgGX = 3, IgGY = 5, L0 = 1e-8, f = 10, murine = true)
    draw(SVG("figure_MITP.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end
