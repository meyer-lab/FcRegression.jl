function loadMixData()
    df = CSV.File(joinpath(dataDir, "lux-mixture.csv"), comment = "#") |> DataFrame!
    df = stack(df, 7:size(df)[2])
    df = dropmissing(df)
    rename!(df, "variable" => "Experiment")
    rename!(df, "value" => "Value")
    df[!, "Value"] = convert.(Float64, df[!, "Value"])

    df[!, "%_1"] ./= 100.0
    df[!, "%_2"] ./= 100.0

    replace!(df."Cell", "hFcgRIIA-131His" => "FcgRIIA-131H")
    replace!(df."Cell", "hFcgRIIB" => "FcgRIIB-232I")
    replace!(df."Cell", "hFcgRIIIA-131Val" => "FcgRIIIA-158V")

    return df
end

function mixNormalExpBatch()
    """ Normalize data without knowing predictions, only by experiment"""
    df = loadMixData()
    meanval = combine(groupby(df, "Experiment"), "Value" => geocmean)
    df = innerjoin(df, meanval, on = "Experiment")
    df[!, "Adjusted"] .= df[!, "Value"] ./ df[!, "Value_geocmean"] .* geocmean(df."Value")
    median(x) = quantile(x, 0.5)
    lower(x) = quantile(x, 0.2)
    upper(x) = quantile(x, 0.8)
    df = combine(groupby(df, ["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]), "Adjusted" => median, "Adjusted" => lower, "Adjusted" => upper)
    rename!(df, "Adjusted_median" => "Value")
    rename!(df, "Adjusted_lower" => "ymin")
    rename!(df, "Adjusted_upper" => "ymax")
    return df
end

function plotMixOriginalData()
    df = mixNormalExpBatch()
    cells = unique(df."Cell")
    pairs = unique(df[!, ["subclass_1", "subclass_2"]])
    lcells = length(cells)
    lpairs = size(pairs, 1)
    pls = Vector(undef, lcells * lpairs)
    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]

    for (i, pairrow) in enumerate(eachrow(pairs))
        for (j, cell) in enumerate(cells)
            IgGXname, IgGYname = pairrow."subclass_1", pairrow."subclass_2"
            ndf = df[(df."Cell" .== cell) .& (df."subclass_1" .== IgGXname) .& (df."subclass_2" .== IgGYname), :]
            pl = plot(
                ndf,
                x = "%_1", 
                y = "Value",
                ymin = "ymin",
                ymax = "ymax",
                color = "Valency", 
                Geom.point,
                Geom.line,
                Geom.errorbar,
                Scale.x_continuous(labels = n -> "$IgGXname $(n*100)%\n$IgGYname $(100-n*100)%"),
                Scale.y_continuous,
                Scale.color_discrete_manual(palette[1], palette[2]),
                Guide.xlabel(""),
                Guide.ylabel("RFU", orientation = :vertical),
                Guide.title("$IgGXname-$IgGYname in $cell"),
            )
            pls[(j - 1) * lpairs + (i - 1) + 1] = pl
        end
    end
    return plotGrid((lcells, lpairs), pls)
end


const measuredRecepExp = Dict("FcgRIIA-131H" => 445141, "FcgRIIB-232I" => 31451, "FcgRIIIA-158V" => 657219)  # geometric mean

function predictMix(dfrow::DataFrameRow, IgGXname, IgGYname, IgGX, IgGY; recepExp = measuredRecepExp)
    IgGC = zeros(size(humanIgG))
    IgGC[IgGXname .== humanIgG] .= IgGX
    IgGC[IgGYname .== humanIgG] .= IgGY

    Kav = importKav(; murine = false, retdf = true)
    Kav = Matrix(Kav[!, [dfrow."Cell"]])
    val = "NewValency" in names(dfrow) ? dfrow."NewValency" : dfrow."Valency"
    res = try
        polyfc(1e-9, KxConst, val, [recepExp[dfrow."Cell"]], IgGC, Kav).Lbound
    catch e
        println(val, [recepExp[dfrow."Cell"]], IgGC, Kav)
        rethrow(e)
    end
    return res
end

predictMix(dfrow::DataFrameRow; recepExp = measuredRecepExp) = predictMix(dfrow, dfrow."subclass_1", dfrow."subclass_2", dfrow."%_1", dfrow."%_2"; recepExp = recepExp)

function predictMix(df::DataFrame; recepExp = measuredRecepExp)
    """ will return another df object """
    df = copy(df)
    df[!, "Predict"] .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictMix(df[i, :]; recepExp = recepExp)
    end
    return df
end


function MixtureFitLoss(df, ValConv::Vector, ExpConv::Vector; logscale = false)
    """ Loss function to be optimized given the conversion factors, called only by MixtureFit() """
    df = copy(df)
    dtype = promote_type(eltype(df."Value"), eltype(ValConv), eltype(ExpConv))
    df[!, "Adjusted"] .= convert.(dtype, df[!, "Value"])
    df[df[!, "Adjusted"] .<= 1.0, "Adjusted"] .= 1.0

    @assert all(df[!, "Adjusted"] .>= 1.0)
    @assert all(df[!, "Predict"] .>= 1.0)
    if any(ValConv .<= 0.0) || any(ExpConv .<= 0.0)
        return Inf, df
    end

    loss = 0.0
    for (iv, val) in enumerate(unique(df."Valency"))
        for (ie, exp) in enumerate(unique(df."Experiment"))
            sdf = @view df[(df."Valency" .== val) .& (df."Experiment" .== exp), :]
            sdf."Adjusted" .*= ValConv[iv] * ExpConv[ie]
            if logscale
                loss += norm(log.(sdf."Adjusted") .- log.(sdf."Predict"), 2) / nrow(sdf)
            else
                loss += norm(sdf."Adjusted" .- sdf."Predict", 2) / nrow(sdf)
            end
        end
    end
    return loss, df
end


function MixtureFit(df; logscale = false)
    """ Two-way fitting for valency and experiment (day) """
    if !("Predict" in names(df))
        df = predictMix(df)
    end
    nv, ne = length(unique(df."Valency")), length(unique(df."Experiment"))
    f(p::Vector, q::Vector) = MixtureFitLoss(df, [1.0; p], q; logscale = logscale)[1]
    f(v::Vector) = f(v[1:(nv - 1)], v[nv:end])
    init_v = ones(nv + ne - 1)
    before_loss = MixtureFitLoss(df, [1.0; init_v[1:(nv - 1)]], init_v[nv:end]; logscale = logscale)[1]
    od = OnceDifferentiable(f, init_v; autodiff = :forward)
    res = optimize(od, init_v, BFGS()).minimizer
    p, q = [1.0; res[1:(nv - 1)]], res[nv:end]
    res = MixtureFitLoss(df, p, q; logscale = logscale)
    return Dict("loss" => res[1], "df" => res[2], 
        "ValConv" => Dict([(name, p[i]) for (i, name) in enumerate(unique(df."Valency"))]), 
        "ExpConv" => Dict([(name, q[i]) for (i, name) in enumerate(unique(df."Experiment"))]))
end

function MixtureCellSeparateFit(df; logscale = false)
    """ Split the cells/receptors and fit valency/exp conv-fac by themselves """
    ndf = nothing
    for cell in unique(df."Cell")
        xdf = MixtureFit(df[df."Cell" .== cell, :]; logscale = logscale)["df"]
        if ndf == nothing
            ndf = xdf
        else
            append!(ndf, xdf)
        end
    end
    return ndf
end


function plotMixContinuous(df; logscale = false)
    """ Use for df with only one pair of IgG subclasses and cell line / Fc receptor """
    df = copy(df)
    @assert length(unique(df."Cell")) == 1
    @assert length(unique(df."subclass_1")) == 1
    @assert length(unique(df."subclass_2")) == 1
    IgGXname = unique(df."subclass_1")[1]
    IgGYname = unique(df."subclass_2")[1]

    x = 0:0.01:1
    df4 = df[(df."Valency" .== 4), :]
    preds4 = [predictMix(df4[1, :], IgGXname, IgGYname, i, 1 - i) for i in x]
    df33 = df[(df."Valency" .== 33), :]
    preds33 = [predictMix(df33[1, :], IgGXname, IgGYname, i, 1 - i) for i in x]

    @assert "Adjusted" in names(df)
    df[!, "Valency"] .= Symbol.(df[!, "Valency"])

    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]
    pl = plot(
        layer(x = x, y = preds4, Geom.line, Theme(default_color = palette[1], line_width = 2px)),
        layer(x = x, y = preds33, Geom.line, Theme(default_color = palette[2], line_width = 2px)),
        layer(df, x = "%_1", y = "Adjusted", color = "Valency", shape = "Experiment"),
        Scale.x_continuous(labels = n -> "$IgGXname $(n*100)%\n$IgGYname $(100-n*100)%"),
        (logscale ? Scale.y_log10(minvalue = 1, maxvalue = 1e6) : Scale.y_continuous),
        Scale.color_discrete_manual(palette[1], palette[2]),
        Guide.xlabel(""),
        Guide.ylabel("RFU", orientation = :vertical),
        Guide.title("$IgGXname-$IgGYname in $(df[1, "Cell"])"),
    )
    return pl
end

function makeMixturePairSubPlots(df; logscale = false)
    cells = unique(df."Cell")
    pairs = unique(df[!, ["subclass_1", "subclass_2"]])
    lcells = length(cells)
    lpairs = size(pairs, 1)
    pls = Vector(undef, lcells * lpairs)

    for (i, pairrow) in enumerate(eachrow(pairs))
        for (j, cell) in enumerate(cells)
            ndf = df[(df."Cell" .== cell) .& (df."subclass_1" .== pairrow."subclass_1") .& (df."subclass_2" .== pairrow."subclass_2"), :]
            pls[(j - 1) * lpairs + (i - 1) + 1] = plotMixContinuous(ndf; logscale = logscale)
        end
    end
    return plotGrid((lcells, lpairs), pls)
end

function plotMixtures()
    """ Plot various fitting options """
    setGadflyTheme()
    res1 = MixtureFit(loadMixData(); logscale = false)
    println("Linear Scale: Fitted ValConv = ", res1["ValConv"], ", ExpConv = ", res1["ExpConv"])
    df1 = res1["df"]
    res2 = MixtureFit(loadMixData(); logscale = true)
    println("Log Scale: Fitted ValConv = ", res2["ValConv"], ", ExpConv = ", res2["ExpConv"])
    df2 = res2["df"]

    df5 = MixtureCellSeparateFit(loadMixData(); logscale = false)
    df6 = MixtureCellSeparateFit(loadMixData(); logscale = true)

    draw(SVG("mix_global_fit_linear.svg", 2500px, 1000px), makeMixturePairSubPlots(df1; logscale = false))
    draw(SVG("mix_global_log.svg", 2500px, 1000px), makeMixturePairSubPlots(df2; logscale = true))
    draw(SVG("mix_cell_split_linear.svg", 2500px, 1000px), makeMixturePairSubPlots(df5; logscale = false))
    draw(SVG("mix_cell_split_log.svg", 2500px, 1000px), makeMixturePairSubPlots(df6; logscale = true))
end
