function figureW(dataType::String; L0 = 1e-9, f = 4, murine::Bool, exp_method = true, fit_ActI = true, legend = true)
    res, loo_res, odf = regressionResult(dataType; L0 = L0, f = f, murine = murine, exp_method = exp_method, fit_ActI = fit_ActI)
    return figureW(res, loo_res, odf, dataType; L0 = L0, f = f, murine = murine, legend = legend)
end

function figureW(res::optResult, loo_res::Vector{optResult}, odf::DataFrame, dataType::String; L0 = 1e-9, f = 4, murine::Bool, legend = true)
    if murine
        df = importDepletion(dataType)
        if dataType == "HIV"
            color = "Label"
        elseif dataType == "Bcell"
            color = "Condition"
        else
            color = "Background"
        end
        shape = "Condition"
    else
        df = importHumanized(dataType)
        color = "Genotype"
        shape = (dataType == "ITP") ? "Condition" : "Concentration"
    end
    @assert all(in(names(odf)).([color, shape]))

    p1 = plotActualvFit(odf, dataType, color, shape; legend = legend)
    p2 = plotCellTypeEffects(df, res, loo_res, dataType; legend = legend, L0 = L0, f = f, murine = murine)
    p3 = plotReceptorActivities(res, loo_res, dataType; murine = murine)
    return p1, p2, p3
end

function plotActualvFit(odf, dataType, colorL::Union{Symbol, String}, shapeL::Union{Symbol, String}; legend = true)
    pl = plot(
        odf,
        x = "Y",
        y = "Fitted",
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


function plotActualvPredict(odf, dataType, colorL::Union{Symbol, String}, shapeL::Union{Symbol, String}; legend = true)
    pl = plot(
        odf,
        x = "Y",
        y = "LOOPredict",
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


function plotCellTypeEffects(df, res, loo_res, dataType; legend = true, L0 = 1e-9, f = 4, murine = true)
    Cell_df = wildtypeWeights(res, df; L0 = L0, f = f, murine = murine)
    Cell_loo = vcat([wildtypeWeights(loo, df) for loo in loo_res]...)
    Cell_conf = combine(groupby(Cell_loo, ["Condition", "Component"]), "Weight" => lower => "ymin", "Weight" => upper => "ymax")
    Cell_df = innerjoin(Cell_df, Cell_conf, on = ["Condition", "Component"])

    pl = plot(
        Cell_df,
        x = "Condition",
        y = "Weight",
        ymin = "ymin",
        ymax = "ymax",
        color = "Component",
        Guide.colorkey(pos = [0.65w, -0.15h]),
        Geom.errorbar,
        Stat.dodge(axis = :x),
        Geom.bar(position = :dodge),
        Scale.x_discrete(levels = unique(Cell_df.Condition)),
        Scale.y_continuous(minvalue = 0.0),
        Scale.color_discrete(levels = unique(Cell_df.Component)),
        Guide.title("Predicted cell type weights for $dataType"),
        style(key_position = legend ? :right : :none, stroke_color = c -> "black"),
    )
    return pl
end

function plotReceptorActivities(res, loo_res, dataType; murine = true)
    ActI_conf = hcat([loo.ActI for loo in loo_res]...)
    ActI_low = lower.(eachslice(ActI_conf, dims = 1))
    ActI_hi = upper.(eachslice(ActI_conf, dims = 1))
    ActI_df = DataFrame(Receptor = (murine ? murineFcgR : humanFcgR), Activity = res.ActI, ymin = ActI_low, ymax = ActI_hi)

    pl = plot(
        ActI_df,
        x = "Receptor",
        y = "Activity",
        ymin = "ymin",
        ymax = "ymax",
        Geom.errorbar,
        Stat.dodge(axis = :x),
        Geom.bar(position = :dodge),
        Scale.x_discrete(),
        Scale.y_continuous(minvalue = 0.0),
        Guide.title("Predicted receptor activities for $dataType"),
        style(bar_spacing = 5mm, stroke_color = c -> "black"),
    )
    return pl
end


function L0fSearchHeatmap(df, dataType, vmax, clmin, clmax; murine = true)
    concs = exp10.(range(clmin, stop = clmax, length = clmax - clmin + 1))
    valencies = [2:vmax;]
    minima = zeros(length(concs), length(valencies))
    for (i, L0) in enumerate(concs)
        for (j, v) in enumerate(valencies)
            Xfc, Xdf, Y = modelPred(df; L0 = L0, f = v, murine = murine)
            fit = fitRegression(Xfc, Xdf, Y; murine = murine)
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
        Guide.title("L<sub>0</sub> and f exploration in $(murine ? "murine" : "human") $dataType data"),
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


function plotSynergy(
    L0,
    f;
    murine = true,
    fit = nothing,
    dataType = "",
    Cellidx = nothing,
    Recepidx = nothing,
    Rbound = false,
    quantity = nothing,
    c1q = false,
    neutralization = false,
)
    Receps = murine ? murineFcgR : humanFcgR
    Kav_df = importKav(; murine = murine, IgG2bFucose = murine, c1q = c1q, retdf = true)
    Kav = Matrix{Float64}(Kav_df[!, Receps])
    if fit === nothing
        title = "Not fit"
    else
        title = "$dataType"
    end

    if Recepidx !== nothing # look at only one receptor
        FcExpr = zeros(length(Receps))
        FcExpr[Recepidx] = importRtot(murine = murine)[Recepidx, Cellidx]
        if murine
            Cell = "$(murineCellTypes[Cellidx])"
            Fc = "$(murineFcgR[Recepidx])"
        else
            Cell = "$(humanCellTypes[Cellidx])"
            Fc = "$(humanFcgR[Recepidx])"
        end
        ylabel = "Activity"
        title = "$title $Cell $Fc"
    elseif Cellidx !== nothing # look at only one cell FcExpr
        FcExpr = importRtot(murine = murine)[:, Cellidx]
        if murine
            Cell = "$(murineCellTypes[Cellidx])"
        else
            Cell = "$(humanCellTypes[Cellidx])"
        end
        ylabel = "Activity"
        title = "$title $Cell"
    else
        FcExpr = importRtot(; murine = murine)
        ylabel = "Depletion"
    end
    if Rbound
        ylabel = "Binding"
    end

    M = synergyGrid(L0, f, FcExpr, Kav; murine = murine, fit = fit, Rbound = Rbound, c1q = c1q, neutralization = neutralization)

    h = collect(Iterators.flatten(M))
    if murine
        S = zeros(length(receptorNamesB1))
        S[1:4] = h[2:5]
        S[5:7] = h[8:10]
        S[8:9] = h[14:15]
        S[10] = h[20]
        S = DataFrame([receptorNamesB1[ii] => S[ii] for ii = 1:length(S)])
    else
        S = zeros(length(humanreceptorNamesB1))
        S[1:3] = h[2:4]
        S[4:5] = h[7:8]
        S[6] = h[12]
        S = DataFrame([humanreceptorNamesB1[ii] => S[ii] for ii = 1:length(S)])
    end

    S = stack(S)

    pl = plot(
        S,
        y = :value,
        x = :variable,
        color = :variable,
        Geom.bar(position = :dodge),
        style(key_position = :none),
        Guide.xlabel("Mixture", orientation = :horizontal),
        Guide.ylabel("Predicted $ylabel Synergy", orientation = :vertical),
        Guide.title("Synergy vs Mixture ($title)"),
    )
    return pl
end
