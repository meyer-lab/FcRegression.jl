function figure_Mmelanoma()
    p1, p2, p3, p4, p5, p6 = figureW("melanoma"; IgGX = 3, IgGY = 5, L0 = 1e-9, f = 6, murine = true)
    draw(SVG("figure_Mmelanoma.svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_MITP()
    p1, p2, p3, p4, p5, p6 = figureW("ITP"; IgGX = 3, IgGY = 5, L0 = 1e-8, f = 10, murine = true)
    draw(SVG("figure_MITP.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_Mblood()
    p1, p2, p3, p4, p5, p6 = figureW("blood"; IgGX = 3, IgGY = 5, L0 = 1e-10, f = 4, murine = true)
    draw(SVG("figure_Mblood.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_Mbone()
    p1, p2, p3, p4, p5, p6 = figureW("bone"; IgGX = 1, IgGY = 2, L0 = 1e-10, f = 4, murine = true)
    draw(SVG("figure_Mbone.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_MHIV()
    p1, p2, p3, p4, p5, p6 = figureW("HIV"; IgGX = 3, IgGY = 5, L0 = 1e-8, f = 24, murine = true)
    draw(SVG("figure_MHIV.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_MBcell()
    p1, p2, p3, p4, p5, p6 = figureW("Bcell", true; IgGX = 1, IgGY = 2, L0 = 1e-10, f = 4, murine = true)
    draw(SVG("figure_MBcell.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end


function figure_Hblood()
    p1, p2, p3, p4, p5, p6 = figureW("blood"; L0 = 1e-8, f = 10, murine = false)
    draw(SVG("figure_Hblood.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_Hbone()
    p1, p2, p3, p4, p5, p6 = figureW("bone"; L0 = 1e-10, f = 4, murine = false)

    draw(SVG("figure_Hspleen.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_Hspleen()
    p1, p2, p3, p4, p5, p6 = figureW("spleen"; L0 = 1e-9, f = 4, murine = false)

    draw(SVG("figure_Hspleen.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end

function figure_HITP()
    p1, p2, p3, p4, p5, p6 = figureW("ITP"; L0 = 1e-9, f = 4, murine = false)
    draw(SVG("figure_HITP.svg", 1000px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
end
