function figureW(dataType; L0 = 1e-9, f = 4, IgGX = 2, IgGY = 3, murine = true, legend = true, Cellidx = nothing, Recepidx = nothing, Rbound = false)
    df = importDepletion(dataType)
    res, odf, Cell_df, ActI_df = regressionResult(dataType; L0 = L0, f = f, murine = murine)

    p1 = plotActualvFit(odf, dataType; legend = legend)
    p2 = plotActualvPredict(odf, dataType; legend = legend)
    p3 = plotCellTypeEffects(Cell_df, dataType; legend = legend)
    p4 = plotReceptorActivities(ActI_df, dataType)
    p5 = plotDepletionSynergy(IgGX, IgGY; L0 = L0, f = f, dataType = dataType, fit = res, Cellidx = Cellidx, Recepidx = Recepidx, Rbound = Rbound)
    p6 = plotSynergy(L0, f; dataType = dataType, fit = res, murine = murine, Cellidx = Cellidx, Recepidx = Recepidx, Rbound = Rbound)

    return p1, p2, p3, p4, p5, p6
end


function plotActualvFit(odf, dataType; legend = true)
    pl = plot(
        odf,
        x = "Y",
        y = "Fitted",
        Geom.point,
        color = "Background",
        shape = "Condition",
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


function plotActualvPredict(odf, dataType; legend = true)
    pl = plot(
        odf,
        x = "Y",
        y = "LOOPredict",
        Geom.point,
        color = "Background",
        shape = "Condition",
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


function plotCellTypeEffects(Cell_df, dataType; legend = true)
    pl = plot(
        Cell_df,
        x = "Condition",
        y = "Weight",
        color = "Component",
        Guide.colorkey(pos = [0.65w, -0.15h]),
        Geom.bar(position = :dodge),
        Scale.x_discrete(levels = unique(Cell_df.Condition)),
        Scale.y_continuous(minvalue = 0.0),
        Scale.color_discrete(levels = unique(Cell_df.Component)),
        Guide.title("Predicted cell type weights for $dataType"),
        style(key_position = legend ? :right : :none),
    )
    return pl
end

function plotReceptorActivities(ActI_df, dataType)
    pl = plot(
        ActI_df,
        x = "Receptor",
        y = "Activity",
        Geom.bar(position = :dodge),
        Scale.x_discrete(),
        Scale.y_continuous(minvalue = 0.0),
        Guide.title("Predicted receptor activities for $dataType"),
        style(bar_spacing = 5mm),
    )
    return pl
end


function L0fSearchHeatmap(df, dataType, vmax, clmin, clmax)
    concs = exp10.(range(clmin, stop = clmax, length = clmax - clmin + 1))
    valencies = [2:vmax;]
    minima = zeros(length(concs), length(valencies))
    for (i, L0) in enumerate(concs)
        for (j, v) in enumerate(valencies)
            Xfc, Xdf, Y = modelPred(df; L0 = L0, f = v)
            fit = fitRegression(Xfc, Xdf, Y)
            minima[i, j] = fit.residual
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
        minima,
        Guide.xlabel("Valencies"),
        Guide.ylabel("L<sub>0</sub> Concentrations"),
        Guide.title("L<sub>0</sub> and f exploration in murine $dataType data"),
        Scale.x_discrete(labels = i -> valencies[i]),
        Scale.y_discrete(labels = i -> concs[i]),
        Scale.color_continuous(minvalue = llim, maxvalue = ulim),
    )
    return pl
end

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

function plotSynergy(L0, f; dataType = nothing, fit = nothing, murine = true, Cellidx = nothing, Recepidx = nothing, Rbound = false, quantity = nothing)
    Kav_df = importKav(; murine = true, IgG2bFucose = true, c1q = false, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, murineFcgR])

    if Recepidx != nothing # look at only one receptor
        FcExpr = zeros(length(murineFcgR))
        FcExpr[Recepidx] = importRtot(; murine = true)[Recepidx, Cellidx]
        ylabel = "Activity"
    elseif Cellidx != nothing # look at only one cell FcExpr
        FcExpr = importRtot(; murine = true)[:, Cellidx]
        ylabel = "Activity"
    else
        FcExpr = importRtot(; murine = true)
    end
    if fit == nothing
        title = "ActI not Fit"
    else
        title = "$dataType"
    end
    if Rbound
        title = "$title Rbound"
    end

    M = synergyGrid(L0, f, FcExpr, Kav; fit = fit, Rbound = Rbound, murine = murine)

    h = collect(Iterators.flatten(M))
    S = zeros(length(receptorNamesB1))
    S[1:4] = h[2:5]
    S[5:7] = h[8:10]
    S[8:9] = h[14:15]
    S[10] = h[20]
    S = DataFrame(Tables.table(S', header = receptorNamesB1))


    S = stack(S)

    pl = plot(
        S,
        y = :value,
        x = :variable,
        color = :variable,
        Geom.bar(position = :dodge),
        style(key_position = :none),
        Guide.xlabel("Mixture", orientation = :vertical),
        Guide.xlabel("Synergy", orientation = :horizontal),
        Guide.title("Synergy vs Mixture ($title)"),
    )
    return pl
end
