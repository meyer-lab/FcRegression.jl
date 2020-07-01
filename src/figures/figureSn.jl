""" This file builds the depletion manuscript, Synergy of ncMO. """
function figureS1()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(1; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS1.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

""" This file builds the depletion manuscript, Synergy of cMO. """
function figureS2()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(2; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS2.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

""" This file builds the depletion manuscript, Synergy of NKs. """
function figureS3()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(3; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS3.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

""" This file builds the depletion manuscript, Synergy of Neu. """
function figureS4()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(4; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS4.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end

""" This file builds the depletion manuscript, Synergy of EO. """
function figureS5()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(5; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS5.svg", 1000px, 1400px), plotGrid((4, 3), [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]))
end
