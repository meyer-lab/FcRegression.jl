function figureS(Cellidx; L0 = 1e-9, f = 4, murine = false)
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


function plotDepletionSynergy(IgGXidx::Int64, IgGYidx::Int64; L0 = 1e-9, f = 4, murine = true,
        fit = nothing, Cellidx = nothing, c1q = false, neutralization = false, RecepIdx = nothing)
    Xname = murine ? murineIgGFucose[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgGFucose[IgGYidx] : humanIgG[IgGYidx]
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    Receps = murine ? murineFcgR : humanFcgR
    ActI = murine ? murineActI : humanActI
    nPoints = 100

    if fit != nothing  # use disease model
        FcExpr = importRtot(; murine = murine)
        title = "Depletion"
        D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr, Kav; murine = murine, fit = fit, ActI = ActI, c1q = c1q)
    elseif Cellidx != nothing
        if murine
            ymax = murineResponse[Cellidx]
        else
            ymax = humanResponse[Cellidx]
        end
        if RecepIdx == nothing  # single cell type
            FcExpr = importRtot(murine = murine)[:, Cellidx]
            title = "Activity"
            D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr, Kav; murine = murine, ActI = ActI)
        else  # bind to one receptor
            FcExpr = zeros(length(Receps), 1)
            FcExpr[RecepIdx, 1] = importRtot(murine = murine)[RecepIdx, Cellidx]
            title = "$(murine ? murineFcgR[RecepIdx] : humanFcgR[RecepIdx]) Binding"
            D1, D2, additive, output = calcSynergy(IgGXidx, IgGYidx, L0, f, FcExpr, Kav; murine = murine)
        end
    else
        @error "Not allowed combination of fit/Cellidx/RecepIdx."
    end

    x = range(0.0, 1.0; length = nPoints)
    @assert length(x) == length(D1)
    @assert length(x) == length(D2)
    @assert length(x) == length(additive)
    @assert length(x) == length(output)

    pl = plot(
        layer(x = x, y = D1, Geom.line, Theme(default_color = colorant"blue", line_width = 1px)),
        layer(x = x, y = D2, Geom.line, Theme(default_color = colorant"orange", line_width = 1px)),
        layer(x = x, y = output, Geom.line, Theme(default_color = colorant"green", line_width = 2px)),
        layer(x = x, y = additive, Geom.line, Theme(default_color = colorant"red", line_width = 3px)),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Guide.xticks(orientation = :horizontal),
        Guide.ylabel("Predicted $title", orientation = :vertical),
        Coord.cartesian(ymin = 0, ymax = ymax),
        Guide.manual_color_key("", ["Predicted", "Additive", "$Xname only", "$Yname only"], ["green", "red", "blue", "orange"]),
        Guide.title("Total predicted effects vs $Xname-$Yname Composition"),
        style(key_position = :inside),
    )
    return pl
end
