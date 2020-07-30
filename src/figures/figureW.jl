function figureW(dataType, intercept = false, preset = false; L0 = 1e-9, f = 4, murine::Bool, ActI = nothing, IgGX = 2, IgGY = 3, legend = true)
    preset_W = nothing
    if murine
        df = importDepletion(dataType)
        color = (dataType == "HIV") ? :Label : :Background
        shape = :Condition
    else
        if preset
            @assert dataType in ["blood", "bone", "ITP"]
            preset_W = fitRegression(importDepletion(dataType), intercept; L0 = L0, f = f, murine = true, ActI = ActI)
            preset_W = preset_W.x
        end
        df = importHumanized(dataType)
        color = :Genotype
        shape = (dataType == "ITP") ? :Condition : :Concentration
    end

    fit, odf, wdf = CVResults(df, intercept, preset_W; L0 = L0, f = f, murine = murine, ActI = ActI)
    @assert all(in(propertynames(odf)).([color, shape]))
    p1 = plotActualvFit(odf, dataType, color, shape; legend = legend)
    p2 = plotActualvPredict(odf, dataType, color, shape; legend = legend)
    p3 = plotCellTypeEffects(wdf, dataType; legend = legend)
    p4 = plotDepletionSynergy(
        IgGX,
        IgGY;
        L0 = L0,
        f = f,
        murine = murine,
        fit = fit,
        c1q = (:C1q in wdf.Component),
        neutralization = (:Neutralization in wdf.Component),
    )
    p5 = L0fSearchHeatmap(df, dataType, 24, -12, -6, murine = murine, ActI = ActI)
    p6 = plotSynergy(L0, f; murine = murine, fit = fit, c1q = (:C1q in wdf.Component), neutralization = (:Neutralization in wdf.Component))

    return p1, p2, p3, p4, p5, p6
end


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
        Guide.title("Actual vs fitted effect for $dataType"),
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
        Guide.title("Actual vs LOO prediction for $dataType"),
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


function L0fSearchHeatmap(df, dataType, vmax, clmin, clmax; murine = true, ActI = nothing)
    concs = exp10.(range(clmin, stop = clmax, length = clmax - clmin + 1))
    valencies = [2:vmax;]
    minimums = zeros(length(concs), length(valencies))
    for (i, L0) in enumerate(concs)
        for (j, v) in enumerate(valencies)
            fit = fitRegression(df; L0 = L0, f = v, murine = murine, ActI = ActI)
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
        Guide.ylabel("L<sub>0</sub> Concentrations"),
        Guide.title("L<sub>0</sub> and f exploration in $(murine ? "murine" : "human") $dataType data"),
        Scale.x_discrete(labels = i -> valencies[i]),
        Scale.y_discrete(labels = i -> concs[i]),
        Scale.color_continuous(minvalue = llim, maxvalue = ulim),
    )
    return pl
end

function plotSynergy(L0, f; murine::Bool, fit = nothing, Cellidx = nothing, quantity = nothing, c1q = false, neutralization = false)
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murine ? murineFcgR : humanFcgR])
    ActI = murine ? murineActI : humanActI

    if Cellidx == nothing #Not using single cell
        FcExpr = importRtot(; murine = murine)
    else #Using single cell
        FcExpr = importRtot(murine = murine)[:, Cellidx]
    end

    M = synergyGrid(L0, f, FcExpr, Kav; murine = murine, fit = fit, ActI = ActI, c1q = c1q)

    h = collect(Iterators.flatten(M))
    if murine
        S = zeros(length(receptorNamesB1))
        S[1:4] = h[2:5]
        S[5:7] = h[8:10]
        S[8:9] = h[14:15]
        S[10] = h[20]
        S = convert(DataFrame, S')
        rename!(S, receptorNamesB1)
    else
        S = zeros(length(humanreceptorNamesB1))
        S[1:3] = h[2:4]
        S[4:5] = h[7:8]
        S[6] = h[12]
        S = convert(DataFrame, S')
        rename!(S, humanreceptorNamesB1)
    end

    S = stack(S)

    pl = plot(S, y = :value, x = :variable, color = :variable, Geom.bar(position = :dodge), style(key_position = :none), Guide.title("Synergy"))
    return pl
end
