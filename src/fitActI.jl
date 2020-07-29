using Optim

function fitActI(dataType; L0 = 1e-10, f = 4, murine::Bool = true, lower::Array{Float64}=nothing, upper::Array{Float64}=nothing)
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
    if murine
        return figureW(dataType; L0 = L0, f = f, murine = true, IgGX = 3, IgGY = 5, ActI = opt.minimizer), opt.minimizer
    else
        return figureW(dataType; L0 = L0, f = f, murine = false, IgGX = 2, IgGY = 3, ActI = opt.minimizer)
    end
end


function figure_fitActI()
    dataTypes = ["ITP", "blood", "bone", "melanoma", "HIV", "Bcell"]
    L0s = [1e-9, 1e-10, 1e-10, 1e-9, 1e-9, 1e-10]
    fs = [6, 5, 4, 6, 4, 4]
    for i in 1:6
        dataType = dataTypes[i]
        pls, w = fitActI(dataType; L0 = L0s[i], f = fs[i], murine = true, lower = [0.0, -5, 0, 0], upper = [5.0, 0, 5, 5])
        p1, p2, p3, p4, p5, p6 = pls
        draw(SVG("figure_fitActI1_M" * dataType * ".svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6], join(string.(w), " "); widths = [5 5 4; 4 6 4]))
        pls, w  = fitActI(dataType; L0 = L0s[i], f = fs[i], murine = true, lower = [0.0, -5, 0, 0], upper = [5.0, 5, 5, 5])
        p1, p2, p3, p4, p5, p6 = pls
        draw(SVG("figure_fitActI2_M" * dataType * ".svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6], join(string.(w), " "); widths = [5 5 4; 4 6 4]))
        pls, w  = fitActI(dataType; L0 = L0s[i], f = fs[i], murine = true, lower = [-5.0, -5, -5, -5], upper = [5.0, 5, 5, 5])
        p1, p2, p3, p4, p5, p6 = pls
        draw(SVG("figure_fitActI3_M" * dataType * ".svg", 1200px, 800px), plotGrid((2, 3), [p1, p2, p3, p4, p5, p6], join(string.(w), " "); widths = [5 5 4; 4 6 4]))
    end
end
