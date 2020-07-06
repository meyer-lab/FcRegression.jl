function figure_Syn_ncMO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(1; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_ncMO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_Syn_cMO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(2; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_cMO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_Syn_NKs()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(3; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_NKs.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_Syn_Neu()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(4; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_Neu.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_Syn_EO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(5; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_EO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end
