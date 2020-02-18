"""Plot figure B5, CD20 blood model"""
function figureB5()
    p1 = plotActualvFit("blood")
    p2 = plotActualvPredict("blood")
    p = plot(p1, p2, layout = (1, 2), size = (1200, 600), dpi = 300)
    savefig(p, "figureB5.pdf")
end
