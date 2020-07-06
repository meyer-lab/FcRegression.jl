""" This file builds the depletion manuscript, Figure 1. """

const receptorNamesB1 =
    Symbol.([
        "IgG1/2a",
        "IgG1/2b",
        "IgG1/3",
        "IgG1/2b-Fucose",
        "IgG2a/2b",
        "IgG2a/3",
        "IgG2a/2b-Fucose",
        "IgG2b/3",
        "IgG2b/2b-Fucose",
        "IgG3/2b-Fucose",
    ])

const humanreceptorNamesB1 = Symbol.(["IgG1/2", "IgG1/3", "IgG1/4", "IgG2/3", "IgG2/4", "IgG3/4"])


"""Figure shows the affect of increasing L0 on binding synergies for each IgG combination"""
function PlotSynGraph(f; murine::Bool, fit = nothing, Cellidx = 2, quantity = nothing, c1q = false, neutralization = false)
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    #ActI = nothing #binding only
    ActI = murine ? murineActI : humanActI
    
    if Cellidx == nothing #Not using single cell
        FcExpr = importRtot(; murine = murine)
    else #Using single cell
        Cell = cellTypes[Cellidx]
        FcExpr = importRtot(murine = murine)[:, Cellidx]
    end
    
    L0 = exp10.(range(-12, stop = -6, length = 20))
    
    if murine
        index = length(receptorNamesB1)
        S = zeros(length(L0), index)
        for (ii, value) in enumerate(L0)
            M = synergyGrid(value, f, FcExpr, Kav; fit = fit, ActI = ActI, c1q = c1q)
            display(M)
            h = collect(Iterators.flatten(M))
            S[ii, 1:4] = h[2:5]
            S[ii, 5:7] = h[8:10]
            S[ii, 8:9] = h[14:15]
            S[ii, 10] = h[20]
        end
        S = convert(DataFrame, S)
        rename!(S, receptorNamesB1)
    else
        index = length(humanreceptorNamesB1)
        S = zeros(length(L0), index)
        for (ii, value) in enumerate(L0)
            M = synergyGrid(value, f, FcExpr, Kav; fit = fit, ActI = ActI, c1q = c1q)
            h = collect(Iterators.flatten(M))
            S[ii, 1:3] = h[2:4]
            S[ii, 4:5] = h[7:8]
            S[ii, 6] = h[12]
        end
        S = convert(DataFrame, S)
        rename!(S, humanreceptorNamesB1)
    end
    S = stack(S)

    pl = plot(
        S,
        x = repeat(L0; outer = [index]),
        y = :value,
        Geom.line,
        color = :variable,
        Guide.colorkey(title = "IgG Combination"),
        Scale.x_log10,
        Guide.xlabel("L0"),
        Guide.ylabel("Synergy"),
        Guide.title("L0 vs Synergy for $Cell"),
        Guide.xticks(),
    )
    return pl
end

""" Figure shows how immune complex valency affects synergy """
function PlotSynValency(L0; murine::Bool, fit = nothing, Cellidx = 2, quantity = nothing, c1q = false, neutralization = false)
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    #ActI = nothing #binding only
    ActI = murine ? murineActI : humanActI
    
    if Cellidx == nothing #Not using single cell
        FcExpr = importRtot(; murine = murine)
    else #Using single cell
        Cell = cellTypes[Cellidx]
        FcExpr = importRtot(murine = murine)[:, Cellidx]
    end
    
    f = range(1, stop = 24)
    
    if murine
        index = length(receptorNamesB1)
        S = zeros(length(f), index)
        for (ii, value) in enumerate(f)
            M = synergyGrid(L0, value, FcExpr, Kav; fit = fit, ActI = ActI, c1q = c1q)
            h = collect(Iterators.flatten(M))
            S[ii, 1:4] = h[2:5]
            S[ii, 5:7] = h[8:10]
            S[ii, 8:9] = h[14:15]
            S[ii, 10] = h[20]
        end
        S = convert(DataFrame, S)
        rename!(S, receptorNamesB1)
    else
        index = length(humanreceptorNamesB1)
        S = zeros(length(f), index)
        for (ii, value) in enumerate(f)
            M = synergyGrid(L0, value, FcExpr, Kav; fit = fit, ActI = ActI, c1q = c1q)
            h = collect(Iterators.flatten(M))
            S[ii, 1:3] = h[2:4]
            S[ii, 4:5] = h[7:8]
            S[ii, 6] = h[12]
        end
        S = convert(DataFrame, S)
        rename!(S, humanreceptorNamesB1)
    end
    S = stack(S)

    pl = plot(
        S,
        x = repeat(f; outer = [index]),
        y = :value,
        Geom.point,
        color = :variable,
        Guide.colorkey(title = "IgG Combination"),
        Guide.xlabel("IC Valency"),
        Guide.ylabel("Synergy"),
        Guide.title("IC Valency vs Synergy for $Cell"),
    )
    return pl
end

function PlotSynvFcrExpr(f, L0; murine::Bool, fit = nothing, Cellidx = 2, quantity = nothing, c1q = false, neutralization = false)
    """ Figure shows how Fc receptor expression affects synergy """
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    #ActI = nothing #binding only
    ActI = murine ? murineActI : humanActI
    
    if Cellidx == nothing #Not using single cell
        FcExpr = importRtot(; murine = murine)
    else #Using single cell
        Cell = cellTypes[Cellidx]
        FcExpr = importRtot(murine = murine)[:, Cellidx]
    end
    
    multiplier = exp10.(range(-2, stop = 0, length = 50))
    S = zeros((50, 10))
    
    if murine
        index = length(receptorNamesB1)
        S = zeros(length(multiplier), index)
        for (ii, value) in enumerate(multiplier)
            M = synergyGrid(L0, f, (FcExpr*value), Kav; fit = fit, ActI = ActI, c1q = c1q)
            h = collect(Iterators.flatten(M))
            S[ii, 1:4] = h[2:5]
            S[ii, 5:7] = h[8:10]
            S[ii, 8:9] = h[14:15]
            S[ii, 10] = h[20]
        end
        S = convert(DataFrame, S)
        rename!(S, receptorNamesB1)
    else
        index = length(humanreceptorNamesB1)
        S = zeros(length(multiplier), index)
        for (ii, value) in enumerate(f)
            M = synergyGrid(L0, f, (FcExpr*value), Kav; fit = fit, ActI = ActI, c1q = c1q)
            h = collect(Iterators.flatten(M))
            S[ii, 1:3] = h[2:4]
            S[ii, 4:5] = h[7:8]
            S[ii, 6] = h[12]
        end
        S = convert(DataFrame, S)
        rename!(S, humanreceptorNamesB1)
    end
    S = stack(S)

    pl = plot(
        S,
        x = repeat(multiplier; outer = [index]),
        y = :value,
        Geom.line,
        Scale.x_log10,
        color = :variable,
        Guide.colorkey(title = "IgG Combination"),
        Guide.xlabel("Fc Expression"),
        Guide.ylabel("Synergy"),
        Guide.title("Fc Expression vs Synergy for $Cell"),
    )
    return pl
end


function figureS(Cellidx; L0 = 1e-9, f = 4, murine = true)
    setGadflyTheme()
    p1 = plotDepletionSynergy(1, 2; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p2 = plotDepletionSynergy(1, 3; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p3 = plotDepletionSynergy(1, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p4 = plotDepletionSynergy(1, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p5 = plotDepletionSynergy(2, 3; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p6 = plotDepletionSynergy(2, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p7 = plotDepletionSynergy(2, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p8 = plotDepletionSynergy(3, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p9 = plotDepletionSynergy(3, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p10 = plotDepletionSynergy(4, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)

    return p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
end
