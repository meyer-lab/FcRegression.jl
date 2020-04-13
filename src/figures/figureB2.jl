""" This file builds the depletion manuscript, Figure 2 (melanoma). """
function figureB2()
    p1, p2, p3, p4 = figureW("melanoma"; L0 = 1e-9, f = 4)

    draw(SVG("figureB2.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end

""" This file builds the depletion manuscript, Figure 3 (ITP). """
function figureB3()
    p1, p2, p3, p4 = figureW("ITP"; L0 = 1e-9, f = 4)

    draw(SVG("figureB3.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end

""" This file builds the depletion manuscript, Figure 4 (blood). """
function figureB4()
    p1, p2, p3, p4 = figureW("blood"; IgGX = 1, IgGY = 3, L0 = 1e-9, f = 4)

    draw(SVG("figureB4.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end

""" This file builds the depletion manuscript, Figure 5 (bone). """
function figureB5()
    p1, p2, p3, p4 = figureW("bone"; IgGX = 1, IgGY = 3, L0 = 1e-9, f = 4)

    draw(SVG("figureB5.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end

""" This file builds the depletion manuscript, Figure 6 (HIV). """
function figureB6()
    p1, p2, p3, p4 = figureW("HIV"; IgGX = 1, IgGY = 2, L0 = 1e-9, f = 4)

    draw(SVG("figureB6.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end
