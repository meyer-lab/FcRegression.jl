function plotActualvFit(odf, dataType, colorL::Symbol, shapeL::Symbol; legend = true)
    pl = plot(
        odf,
        x = :Y,
        y = :Fitted,
        Geom.point,
        color = colorL,
        shape = shapeL,
        Guide.colorkey(),
        Guide.shapekey(),
        Scale.y_continuous(minvalue = 0.0, maxvalue = 1.0),
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("Fitted effect"),
        Guide.title("Actual effect vs fitted effect for $dataType"),
        style(point_size = 5px, key_position = legend ? :right : :none),
    )
    return pl
end


function plotActualvPredict(odf, dataType, colorL::Symbol, shapeL::Symbol; legend = true)
    pl = plot(
        odf,
        x = :Y,
        y = :LOOPredict,
        Geom.point,
        color = colorL,
        shape = shapeL,
        Guide.colorkey(),
        Guide.shapekey(),
        Geom.abline(color = "red"),
        Guide.xlabel("Actual effect"),
        Guide.ylabel("LOO predicted effect"),
        Guide.title("Actual effect vs LOO prediction for $dataType"),
        style(point_size = 5px, key_position = legend ? :right : :none),
    )
    return pl
end


function plotCellTypeEffects(wdf, dataType; legend = true)
    pl = plot(
        wdf,
        x = :Condition,
        y = :Median,
        color = :Component,
        Guide.colorkey(pos = [0.05w, -0.3h]),
        Geom.bar(position = :dodge),
        Scale.x_discrete(levels = unique(wdf.Condition)),
        Scale.y_continuous(minvalue = 0.0),
        Scale.color_discrete(levels = unique(wdf.Component)),
        ymin = :Q10,
        ymax = :Q90,
        Geom.errorbar,
        Stat.dodge,
        Guide.title("Predicted cell type weights for $dataType"),
        style(key_position = legend ? :right : :none),
    )
    return pl
end


function plotDepletionSynergy(IgGXidx::Int64, IgGYidx::Int64, fit::fitResult; L0, f, murine::Bool, c1q = false, neutralization = false)
    Xname = murine ? murineIgG[IgGXidx] : humanIgG[IgGXidx]
    Yname = murine ? murineIgG[IgGYidx] : humanIgG[IgGYidx]
    Kav_df = importKav(; murine = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    FcExpr = importRtot(; murine = murine)
    ActI = murine ? murineActI : humanActI

    nPoints = 100
    IgGC = zeros(Float64, size(Kav, 1), nPoints)

    IgGC[IgGYidx, :] = range(eps(), eps(); length = nPoints)
    IgGC[IgGXidx, :] = range(1.0, 1.0; length = nPoints)
    X1 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI, Mix = false)  # size: celltype * nPoints

    IgGC[IgGXidx, :] = range(eps(), eps(); length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 1.0; length = nPoints)
    X2 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI, Mix = false)  # size: celltype * nPoints

    IgGC[IgGXidx, :] = range(0.0, 1.0; length = nPoints)
    IgGC[IgGYidx, :] = range(1.0, 0.0; length = nPoints)
    X = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI)  # size: celltype * nPoints

    if c1q
        X = vcat(X, Kav_df[!, :C1q]' * IgGC)
        X1 = vcat(X1, Kav_df[!, :C1q]' * IgGC)
        X2 = vcat(X2, Kav_df[!, :C1q]' * IgGC)
    end
    @assert size(X, 1) == length(fit.x)
    @assert size(X1, 1) == length(fit.x)
    @assert size(X2, 1) == length(fit.x)

    output = exponential(Matrix(X'), fit)
    D1 = exponential(Matrix(X1'), fit)
    D2 = reverse(exponential(Matrix(X2'), fit))

    pl = plot(
        layer(x = IgGC[IgGXidx, :], y = output, Geom.line, Theme(default_color = colorant"green")),
        layer(x = IgGC[IgGXidx, :], y = D1 + D2, Geom.line, Theme(default_color = colorant"red")),
        layer(x = IgGC[IgGXidx, :], y = D1, Geom.line, Theme(default_color = colorant"blue")),
        layer(x = IgGC[IgGXidx, :], y = D2, Geom.line, Theme(default_color = colorant"yellow")),
        Scale.x_continuous(labels = n -> "$Xname $(n*100)%\n$Yname $(100-n*100)%"),
        Guide.ylabel("Predicted Depletion"),
        Guide.manual_color_key("", ["Predicted", "Additive", "$Xname only", "$Yname only"], ["green", "red", "blue", "yellow"]),
        Guide.title("Total predicted effects vs $Xname-$Yname Composition"),
        style(key_position = :inside),
    )
    return pl
end


function createHeatmap(df, dataType, vmax, clmin, clmax; murine = true)
    concs = exp10.(range(clmin, stop = clmax, length = clmax - clmin + 1))
    valencies = [2:vmax;]
    minimums = zeros(length(concs), length(valencies))
    for (i, L0) in enumerate(concs)
        for (j, v) in enumerate(valencies)
            fit = fitRegression(df; L0 = L0, f = v, murine = murine)
            minimums[i, j] = fit.r
        end
    end
    if dataType == "HIV"
        llim = 0
        ulim = 0.1
    else
        llim = nothing
        ulim = nothing
    end
    pl = spy(
        minimums,
        Guide.xlabel("Valencies"),
        Guide.ylabel("L0 Concentrations"),
        Guide.title("L_0 and f exploration in $(murine ? "murine" : "human") $dataType data"),
        Scale.x_discrete(labels = i -> valencies[i]),
        Scale.y_discrete(labels = i -> concs[i]),
        Scale.color_continuous(minvalue = llim, maxvalue = ulim),
    )
    return pl
end

function plotSynergy(fit::fitResult; L0, f, murine::Bool, c1q = false, neutralization = false)
    Kav_df = importKav(; murine = murine, IgG2bFucose = true, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    FcExpr = importRtot(; murine = murine)
    ActI = murine ? murineActI : humanActI

    nPoints = 100
    M = zeros(size(Kav)[1], size(Kav)[1])

    for i = 1:size(Kav)[1]
        for j = 1:(i - 1)
            IgGC = zeros(Float64, size(Kav, 1), nPoints)
            IgGC[j, :] .= eps()
            IgGC[i, :] .= 1
            X1 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI, Mix = false)  # size: celltype * nPoints

            IgGC[i, :] .= eps()
            IgGC[j, :] .= 1
            X2 = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI, Mix = false)  # size: celltype * nPoints

            IgGC[i, :] = range(0.0, 1.0; length = nPoints)
            IgGC[j, :] = range(1.0, 0.0; length = nPoints)
            X = polyfc_ActV(L0, KxConst, f, FcExpr, IgGC, Kav, ActI)  # size: celltype * nPoints
            if c1q
                X = vcat(X, Kav_df[!, :C1q]' * IgGC)
                X1 = vcat(X1, Kav_df[!, :C1q]' * IgGC)
                X2 = vcat(X2, Kav_df[!, :C1q]' * IgGC)
            end
            @assert size(X, 1) == length(fit.x)
            @assert size(X1, 1) == length(fit.x)
            @assert size(X2, 1) == length(fit.x)
            output = exponential(Matrix(X'), fit)
            D1 = exponential(Matrix(X1'), fit)
            D2 = reverse(exponential(Matrix(X2'), fit))
            additive = D1 + D2
            synergy = sum((output - additive) / nPoints)
            M[i, j] = synergy
        end
    end
    
    h = collect(Iterators.flatten(M))
    if murine
        S = zeros(10)
        S[1:4] = h[2:5]
        S[5:7] = h[8:10]
        S[8:9] = h[14:15]
        S[10] = h[20]
        S = convert(DataFrame, S')
        rename!(S, Symbol.(receptorNamesB1()))
    else
        S = zeros(6)
        S[1:3] = h[2:4]
        S[4:5] = h[7:8]
        S[6] = h[12]
        S = convert(DataFrame, S')
        rename!(S, Symbol.(humanreceptorNamesB1()))
    end

    S = stack(S)

    pl = plot(S, y = :value, x = :variable, color = :variable, Geom.bar(position = :dodge), Theme(key_position = :none), Guide.title("Synergy"))
    return pl
end

function figureW(dataType, intercept = false, preset = false; L0 = 1e-9, f = 4, murine::Bool, IgGX = 2, IgGY = 3, legend = true)
    preset_W = nothing
    if murine
        df = importDepletion(dataType)
        color = (dataType == "HIV") ? :Label : :Background
        shape = :Condition
    else
        if preset
            @assert dataType in ["blood", "bone", "ITP"]
            preset_W = fitRegression(importDepletion(dataType), intercept; L0 = L0, f = f, murine = true)
            preset_W = preset_W.x
        end
        df = importHumanized(dataType)
        color = :Genotype
        shape = (dataType == "ITP") ? :Condition : :Concentration
    end

    fit, odf, wdf = CVResults(df, intercept, preset_W; L0 = L0, f = f, murine = murine)
    @assert all(in(propertynames(odf)).([color, shape]))
    p1 = plotActualvFit(odf, dataType, color, shape; legend = legend)
    p2 = plotActualvPredict(odf, dataType, color, shape; legend = legend)
    p3 = plotCellTypeEffects(wdf, dataType; legend = legend)
    p4 = plotDepletionSynergy(
        IgGX,
        IgGY,
        fit;
        L0 = L0,
        f = f,
        murine = murine,
        c1q = (:C1q in wdf.Component),
        neutralization = (:Neutralization in wdf.Component),
    )
    p5 = createHeatmap(df, dataType, 24, -12, -6, murine = murine)
    p6 = plotSynergy(fit; L0 = L0, f = f, murine = murine, c1q = (:C1q in wdf.Component), neutralization = (:Neutralization in wdf.Component))

    return p1, p2, p3, p4, p5, p6
end
