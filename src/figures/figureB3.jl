""" This file builds the depletion manuscript, Figure 3 (melanoma). """

function figureB3()
    p1, p2, p3, p4 = figureW("melonoma"; L0 = 1e-9, f = 4)

    draw(PDF("figureB3.pdf", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end
