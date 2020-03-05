function plotActualvFit(odf, dataType)
    fit_res = lm(@formula(Fitted ~ Y), odf)
    intercept, slope = coef(fit_res)
    r_sq = r2(fit_res)

    pl = plot(
        odf,
        x = :Y,
        y = :Fitted,
        Geom.point,
        color = :Condition,
        shape = :Background,
        Guide.colorkey(pos = [0.05w, 0.28h]),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        slope = [slope],
        intercept = [intercept],
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("Fitted effect"),
        Guide.title("Actual effect vs fitted effect for $dataType"),
        Theme(point_size = 5px),
    )
    return pl
end


function plotActualvPredict(odf, dataType)
    fit_res = lm(@formula(Fitted ~ Y), odf)
    intercept, slope = coef(fit_res)
    r_sq = r2(fit_res)

    pl = plot(
        odf,
        x = :Y,
        y = :LOOPredict,
        Geom.point,
        color = :Condition,
        shape = :Background,
        Guide.colorkey(pos = [0.05w, 0.5h]),
        slope = [slope],
        intercept = [intercept],
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("LOO predicted effect"),
        Guide.title("Actual effect vs LOO predicted for $dataType"),
    )
    return pl
end


function plotCellTypeEffects(wdf, dataType)
    wdf.ymin = wdf.Weight .- wdf.BtpStdev
    wdf.ymax = wdf.Weight .+ wdf.BtpStdev

    pl = plot(
        wdf,
        x = :Condition,
        y = :Weight,
        color = :Component,
        Guide.colorkey(pos = [0.05w, -0.28h]),
        Geom.bar(position = :dodge),
        Scale.x_discrete(levels = unique(wdf.Condition)),
        Scale.color_discrete(levels = unique(wdf.Component)),
        ymin = :ymin,
        ymax = :ymax,
        Geom.errorbar,
        Stat.dodge,
        Guide.title("Weights of IgG and Cell Type in wt for $dataType"),
    )
    return pl
end


function plotDepletionSynergy(IgGXidx::Int64, IgGYidx::Int64, weights::Vector, L0, f; murine = true, nPoints = 40, c1q = false)
    Xname = murine ? murineIgG[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgG[IgGYidx] : humanIgG[IgGYidx]
    Kav_df = importKav(; murine = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    FcExpr = importRtot(; murine = murine)
    ActI = murine ? murineActI : humanActI

    IgGC = zeros(Float64, size(Kav, 1), nPoints)
    IgGC[IgGXidx, :] = range(0.0, 1.0; length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 0.0; length = nPoints)
    X = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI)  # size: celltype * nPoints
    if c1q
        X = vcat(X, Kav_df[!, :C1q]' * IgGC)
    end
    @assert size(X, 1) == length(weights)
    output = exponential(Matrix(X'), weights)

    pl = plot(
        layer(x = IgGC[IgGXidx, :], y = output, Geom.line, Theme(default_color = colorant"green")),
        layer(x = [0, 1], y = [output[1], output[end]], Geom.line, Theme(default_color = colorant"red")),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Guide.manual_color_key("", ["Predicted", "Linear Addition"], ["green", "red"]),
        Guide.title("Total predicted effects vs $Xname-$Yname Composition"),
        Theme(key_position = :inside),
    )
    return pl
end


function figureW(dataType; L0 = 1e-9, f = 4)
    (fit_w, odf, wdf) = CVResults(dataType; L0 = L0, f = f)
    p1 = plotActualvFit(odf, dataType)
    p2 = plotActualvPredict(odf, dataType)
    p3 = plotCellTypeEffects(wdf, dataType)
    p4 = plotDepletionSynergy(2, 3, fit_w, 1e-9, 4; c1q = (:C1q in unique(wdf.Component)))

    return p1, p2, p3, p4
end
