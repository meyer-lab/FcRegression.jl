using Optim

function fitActI(dataType; L0 = 1e-10, f = 4, murine::Bool = true, lower=nothing, upper=nothing)
    df = murine ? importDepletion(dataType) : importHumanized(dataType)
    func = ActI -> fitRegression(df; L0 = L0, f = f, murine=murine, ActI = ActI).r
    if lower == nothing && upper == nothing
        if murine
            lower = [0.0, -10.0, 0.0, 0.0]
            upper = [10.0, 10.0, 10.0, 10.0]
        else
            lower = [0.0, 0.0, 0.0, -10.0, -10.0, 0.0, 0.0, 0.0, 0.0]
            upper = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
        end
    end
    opt = optimize(func, lower, upper, Float64.(murine ? murineActI : humanActI))
    #display(opt)
    display(opt.minimizer)
    if murine
        return figureW(dataType; L0 = L0, f = f, murine = true, IgGX = 3, IgGY = 5, ActI = opt.minimizer), opt.minimizer
    else
        return figureW(dataType; L0 = L0, f = f, murine = false, IgGX = 2, IgGY = 3, ActI = opt.minimizer)
    end
end


function figure_fitActI()
    for dataType in ["ITP", "blood", "bone", "melanoma", "HIV", "Bcell"]
        pls, w = fitActI(dataType; L0 = 1e-10, f = 6, murine = true, lower = [0.0, -10.0, 0.0, 0.0], upper = [10.0, 0.0, 10.0, 10.0])
        p1, p2, p3, p4, p5, p6 = pls
        draw(SVG("figure_fitActI1_M" * dataType * ".svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6], join(string.(w), " "); widths = [5 5 4; 4 6 4]))
        pls, w  = fitActI(dataType; L0 = 1e-10, f = 6, murine = true, lower = [0.0, -10.0, 0.0, 0.0], upper = [10.0, 10.0, 10.0, 10.0])
        p1, p2, p3, p4, p5, p6 = pls
        draw(SVG("figure_fitActI2_M" * dataType * ".svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6], join(string.(w), " "); widths = [5 5 4; 4 6 4]))
        pls, w  = fitActI(dataType; L0 = 1e-10, f = 6, murine = true, lower = [-10.0, -10.0, -10.0, -10.0], upper = [10.0, 10.0, 10.0, 10.0])
        p1, p2, p3, p4, p5, p6 = pls
        draw(SVG("figure_fitActI3_M" * dataType * ".svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6], join(string.(w), " "); widths = [5 5 4; 4 6 4]))
    end
    #=for dataType in ["blood", "spleen", "bone", "ITP"]
        p1, p2, p3, p4, p5, p6 = fitActI(dataType; L0 = 1e-10, f = 6, murine = false)
        draw(SVG("figure_fitActI_H" * dataType * ".svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6]; widths = [5 5 4; 4 6 4]))
    end=#
end
