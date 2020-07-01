""" This file builds the depletion manuscript, Synergy of ncMO. """
function figure_SynncMO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(1; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_SynncMO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

""" This file builds the depletion manuscript, Synergy of cMO. """
function figure_SyncMO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(2; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_SyncMO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

""" This file builds the depletion manuscript, Synergy of NKs. """
function figure_SynNKs()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(3; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_SynNKs.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

""" This file builds the depletion manuscript, Synergy of Neu. """
function figure_SynNeu()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(4; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_SynNeu.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

""" This file builds the depletion manuscript, Synergy of EO. """
function figure_SynEO()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(5; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figure_SynEO.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end
