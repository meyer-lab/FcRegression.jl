""" This file builds the depletion manuscript, Figure 2 (melanoma). """
function figureB2()
    p1, p2, p3, p4 = figureW("melanoma"; L0 = 1e-9, f = 4)

    draw(SVG("figureB2.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end

""" This file builds the depletion manuscript, Figure 3 (ITP). """
function figureB3()
    p1, p2, p3, p4 = figureW("ITP"; L0 = 1e-9, f = 4)

    draw(PDF("figureB3.pdf", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end

""" This file builds the depletion manuscript, Figure 4 (blood). """
function figureB4()
    p1, p2, p3, p4 = figureW("blood"; L0 = 1e-9, f = 4)

    draw(PDF("figureB4.pdf", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end

""" This file builds the depletion manuscript, Figure 5 (bone). """
function figureB5()
    p1, p2, p3, p4 = figureW("bone"; L0 = 1e-9, f = 4)

    draw(PDF("figureB5.pdf", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end

""" This file builds the depletion manuscript, Figure 6 (HIV). """
function figureB6()
    p1, p2, p3, p4 = figureW("HIV"; L0 = 1e-9, f = 4)

    draw(PDF("figureB5.pdf", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end
