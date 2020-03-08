""" This file builds the depletion manuscript, Figure 2 (ITP). """

function figureB2()
    p1, p2, p3, p4 = figureW("ITP"; L0 = 1e-9, f = 4)

    draw(SVG("figureB2.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end
