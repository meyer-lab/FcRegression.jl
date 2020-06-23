""" This file builds the depletion manuscript, Figure 2 (melanoma). """
function figureB2()
    p1, p2, p3, p4, p5, p6 = figureW("melanoma"; IgGX = 2, IgGY = 3, L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureB2.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Figure 3 (ITP). """
function figureB3()
    p1, p2, p3, p4, p5, p6 = figureW("ITP"; IgGX = 1, IgGY = 2, L0 = 1e-8, f = 10, murine = true)
    draw(SVG("figureB3.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Figure 4 (blood). """
function figureB4()
    p1, p2, p3, p4, p5, p6 = figureW("blood"; IgGX = 1, IgGY = 2, L0 = 1e-10, f = 4, murine = true)
    draw(SVG("figureB4.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Figure 5 (bone). """
function figureB5()
    p1, p2, p3, p4, p5, p6 = figureW("bone"; IgGX = 1, IgGY = 2, L0 = 1e-10, f = 4, murine = true)
    draw(SVG("figureB5.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Figure 6 (HIV). """
function figureB6()
    p1, p2, p3, p4, p5, p6 = figureW("HIV"; IgGX = 1, IgGY = 2, L0 = 1e-8, f = 24, murine = true)
    draw(SVG("figureB6.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" Figure 7 (bone) from Lux C57BL/6 data """
function figureB7()
    p1, p2, p3, p4, p5, p6 = figureW("Bcell", true; IgGX = 2, IgGY = 4, L0 = 1e-10, f = 4, murine = true)
    draw(SVG("figureB7.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 11 (blood). """
function figureB11()
    p1, p2, p3, p4, p5, p6 = figureW("blood"; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figureB11.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 12 (spleen). """
function figureB12()
    p1, p2, p3, p4, p5, p6 = figureW("spleen"; L0 = 1e-9, f = 4, murine = false)

    draw(SVG("figureB12.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 13 (HIV). """
function figureB13()
    p1, p2, p3, p4, p5, p6 = figureW("bone"; L0 = 1e-9, f = 4, murine = false)

    draw(SVG("figureB13.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 14 (blood). """
function figureB14()
    p1, p2, p3, p4, p5, p6 = figureW("blood", true, true; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figureB14.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 15 (bone). """
function figureB15()
    p1, p2, p3, p4, p5, p6 = figureW("bone", true, true; L0 = 1e-9, f = 4, murine = false)

    draw(SVG("figureB15.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 16 (ITP) from Schwab 2015. """
function figureB16()
    p1, p2, p3, p4, p5, p6 = figureW("ITP"; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figureB16.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 17 (bone) from Schwab 2015. """
function figureB17()
    p1, p2, p3, p4, p5, p6 = figureW("ITP", true, true; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figureB17.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Synergy of ncMO. """
function figureS1()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(1; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS1.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p10 p10]))
end

""" This file builds the depletion manuscript, Synergy of cMO. """
function figureS2()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(2; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS2.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p10 p10]))
end

""" This file builds the depletion manuscript, Synergy of NKs. """
function figureS3()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(3; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS3.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p10 p10]))
end

""" This file builds the depletion manuscript, Synergy of Neu. """
function figureS4()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(4; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS4.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p10 p10]))
end

""" This file builds the depletion manuscript, Synergy of EO. """
function figureS5()
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = figureS(5; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureS5.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p10 p10]))
end
