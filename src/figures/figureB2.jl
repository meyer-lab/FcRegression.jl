""" This file builds the depletion manuscript, Figure 2 (melanoma). """
function figureB2()
    p1, p2, p3, p4, p5, p6 = figureW("melanoma"; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureB2.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Figure 3 (ITP). """
function figureB3()
    p1, p2, p3, p4, p5, p6 = figureW("ITP"; L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureB3.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Figure 4 (blood). """
function figureB4()
    p1, p2, p3, p4, p5, p6 = figureW("blood"; IgGX = 1, IgGY = 3, L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureB4.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Figure 5 (bone). """
function figureB5()
    p1, p2, p3, p4, p5, p6 = figureW("bone"; IgGX = 1, IgGY = 3, L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureB5.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, Figure 6 (HIV). """
function figureB6()
    p1, p2, p3, p4, p5, p6 = figureW("HIV"; IgGX = 1, IgGY = 3, L0 = 1e-9, f = 4, murine = true)
    draw(SVG("figureB6.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p5 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 11 (blood). """
function figureB11()
    p1, p2, p3, p4, p5, p6 = figureW("blood"; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figureB11.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p6 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 12 (spleen). """
function figureB12()
    p1, p2, p3, p4, p5, p6 = figureW("spleen"; L0 = 1e-9, f = 4, murine = false)

    draw(SVG("figureB12.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p6 p6]))
end

""" This file builds the depletion manuscript, humanized mice, Figure 13 (HIV). """
function figureB13()
    p1, p2, p3, p4, p5, p6 = figureW("bone marrow"; L0 = 1e-9, f = 4, murine = false)

    draw(SVG("figureB13.svg", 1000px, 800px), gridstack([p1 p2 p3; p4 p6 p6]))
end
