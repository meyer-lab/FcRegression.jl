function figure_MSyn_ncMO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(1; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_MncMO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_MSyn_cMO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(2; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_McMO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_MSyn_NKs()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(3; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_MNKs.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_MSyn_Neu()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(4; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_MNeu.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_MSyn_EO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(5; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_Syn_MEO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_HSyn_ncMO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(1; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figure_Syn_HncMO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_HSyn_cMO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(2; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figure_Syn_HcMO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_HSyn_NKs()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(3; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figure_Syn_HNKs.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_HSyn_Neu()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(4; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figure_Syn_HNeu.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

function figure_HSyn_EO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(5; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figure_Syn_HEO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end
