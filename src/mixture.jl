using Dierckx
using MultivariateStats
using Impute
using StatsBase
using GLM

function loadMixData(fn = "lux_mixture_mar2021.csv"; avg = false)
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame

    #appends average column
    if avg
        av_df = copy(df)
        for col in eachcol(av_df)
            replace!(col,missing => 0)
        end
        av = (sum(eachcol(av_df[:,7:ncol(av_df)])) ./(ncol(av_df)-6))
        df = df[:,1:6]
        df[!, :average] = av
    end

    df = stack(df, 7:size(df)[2])
    df = dropmissing(df)
    rename!(df, "variable" => "Experiment")
    rename!(df, "value" => "Value")
    df[!, "Value"] = convert.(Float64, df[!, "Value"])

    df[!, "%_1"] ./= 100.0
    df[!, "%_2"] ./= 100.0

    replace!(df."Cell", "CHO-hFcgRIIA-131His" => "FcgRIIA-131H")
    replace!(df."Cell", "CHO-hFcgRIIB" => "FcgRIIB-232I")
    replace!(df."Cell", "CHO-hFcgRIIIA-131Val" => "FcgRIIIA-158V")
    replace!(df."Cell", "CHO-FcgRIA" => "FcgRI")
    replace!(df."Cell", "CHO-hFcgRIIA-131Arg" => "FcgRIIA-131R")
    replace!(df."Cell", "CHO-hFcgRIIIA-158Phe" => "FcgRIIIA-158F")

    return sort!(df, ["Valency", "Cell", "subclass_1", "subclass_2", "Experiment", "%_2"])
end

function mixNormalExpBatch(df = loadMixData())
    """ Normalize data without knowing predictions, only by experiment"""
    #meanval = combine(groupby(df, "Experiment"), "Value" => geocmean)
    #df = innerjoin(df, meanval, on = "Experiment")
    df[!, "Adjusted"] .= df[!, "Value"] #./ df[!, "Value_geocmean"] .* geocmean(df."Value")
    median(x) = quantile(x, 0.5)
    lower(x) = quantile(x, 0.2)
    upper(x) = quantile(x, 0.8)
    df = combine(
        groupby(df, ["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]),
        "Adjusted" => median,
        "Adjusted" => lower,
        "Adjusted" => upper,
    )
    rename!(df, "Adjusted_median" => "Value")
    rename!(df, "Adjusted_lower" => "ymin")
    rename!(df, "Adjusted_upper" => "ymax")
    return df
end

function plotMixOriginalData(df = loadMixData())
    df = mixNormalExpBatch(df)
    cells = unique(df."Cell")
    pairs = unique(df[!, ["subclass_1", "subclass_2"]])
    lcells = length(cells)
    lpairs = size(pairs, 1)
    pls = Vector(undef, lcells * lpairs)
    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]
    ymax = Dict(
        "FcgRI" => 8e3,
        "FcgRIIA-131H" => 2.5e4,
        "FcgRIIA-131R" => 2.5e4,
        "FcgRIIB-232I" => 3e3,
        "FcgRIIIA-158F" => 2e4,
        "FcgRIIIA-158V" => 1.5e4,
    )

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
                Scale.y_continuous(; maxvalue = ymax[cell]),
                Scale.color_discrete_manual(palette[1], palette[2]),
                Guide.xlabel("", orientation = :horizontal),
                Guide.ylabel("RFU", orientation = :vertical),
                Guide.xticks(orientation = :horizontal),
                Guide.title("$IgGXname-$IgGYname in $cell"),
            )
            pls[(j - 1) * lpairs + (i - 1) + 1] = pl
        end
    end
    return plotGrid((lcells, lpairs), pls; sublabel = false)
end

function mixEC50()
    df = mixNormalExpBatch()
    cells = unique(df."Cell")
    pairs = unique(df[!, ["subclass_1", "subclass_2"]])
    lcells = length(cells)
    lpairs = size(pairs, 1)

    Bound = Vector(undef, lcells * lpairs)
    Ka = Vector(undef, lcells * lpairs)
    Combos = Vector(undef, lcells * lpairs)
    Combos = Vector(undef, lcells * lpairs)
    Cells = Vector(undef, lcells * lpairs)

    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3]]
    Kav = importKav(; murine = false, retdf = true)

    for (i, pairrow) in enumerate(eachrow(pairs))
        for (j, cell) in enumerate(cells)
            IgGXname, IgGYname = pairrow."subclass_1", pairrow."subclass_2"
            ndf = df[(df."Cell" .== cell) .& (df."subclass_1" .== IgGXname) .& (df."Valency" .== 4) .& (df."subclass_2" .== IgGYname), :]
            sort!(ndf, ["%_1"])
            y = ndf["Value"]
            x = ndf["%_1"]
            sp = Spline1D(x, y)
            x = 0:0.01:1.0
            y = sp(x)
            EC50value = 0.5 * maximum(y)
            diff = y .- EC50value
            Value, EC50index = findmin(abs.(diff))
            if Kav[j, (2^j) + 1] < Kav[j + 1, (2^j) + 1]
                Ka[(j - 1) * lpairs + (i - 1) + 1] = Kav[j + 1, (2^j) + 1]
                Cells[(j - 1) * lpairs + (i - 1) + 1] = string(IgGYname, " ", cells[j])
            else
                Ka[(j - 1) * lpairs + (i - 1) + 1] = Kav[j, (2^j) + 1]
                Cells[(j - 1) * lpairs + (i - 1) + 1] = string(IgGXname, " ", cells[j])
            end
            Bound[(j - 1) * lpairs + (i - 1) + 1] = Value
            Combos[(j - 1) * lpairs + (i - 1) + 1] = "$IgGXname/$IgGYname"
        end
    end

    p1 = plot(
        x = Ka,
        y = Bound,
        color = Combos,
        shape = Cells,
        Geom.point,
        Scale.x_log10,
        Scale.y_continuous,
        Scale.color_discrete_manual(palette[1], palette[2]),
        Guide.xlabel("Kav"),
        Guide.ylabel("EC50 Rbound", orientation = :vertical),
        style(point_size = 5px, key_position = :right),
    )
    return p1
end


const measuredRecepExp = Dict(
    "FcgRI" => 101493.689,
    "FcgRIIA-131H" => 1006302.484,
    "FcgRIIA-131R" => 190432.6753,
    "FcgRIIB-232I" => 75085.07599,
    "FcgRIIIA-158F" => 634324.0675,
    "FcgRIIIA-158V" => 979451.9884,
)  # geometric mean


function R2(Actual, Predicted)
    df = DataFrame(A=log10.(Actual), B=log10.(Predicted))
    ols = lm(@formula(B ~ A + 0), df)
    display(ols)
    R2 = r2(ols)
    return R2
end


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

predictMix(dfrow::DataFrameRow; recepExp = measuredRecepExp) =
    predictMix(dfrow, dfrow."subclass_1", dfrow."subclass_2", dfrow."%_1", dfrow."%_2"; recepExp = recepExp)

function predictMix(df::DataFrame; recepExp = measuredRecepExp)
    """ will return another df object """
    df = copy(df)
    df[!, "Predict"] .= 1.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictMix(df[i, :]; recepExp = recepExp)
    end
    df[df."Predict" .< 1.0, "Predict"] .= 1.0
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
    return Dict(
        "loss" => res[1],
        "df" => res[2],
        "ValConv" => Dict([(name, p[i]) for (i, name) in enumerate(unique(df."Valency"))]),
        "ExpConv" => Dict([(name, q[i]) for (i, name) in enumerate(unique(df."Experiment"))]),
    )
end

function MixtureCellSeparateFit(df; logscale = false)
    """ Split the cells/receptors and fit valency/exp conv-fac by themselves """
    ndf = DataFrame()
    for cell in unique(df."Cell")
        append!(ndf, MixtureFit(df[df."Cell" .== cell, :]; logscale = logscale)["df"])
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

    if !("Adjusted" in names(df))
        df[!, "Adjusted"] .= df[!, "Value"]
    end
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
        Guide.xticks(orientation = :horizontal),
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


function PCAData(; cutoff = 0.9)
    df = loadMixData()
    exps_sets = [["6/23/20", "6/30/20", "7/14/20", "7/23/20", "9/11/20"], ["5/15/20", "5/20/20", "5/28/20", "6/2/20", "9/2/20"]]
    retdf = DataFrame()

    for (i, exps) in enumerate(exps_sets)
        ndf = df[in(exps).(df."Experiment"), :]
        widedf = unstack(ndf, ["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"], "Experiment", "Value")
        widedf = coalesce.(widedf, 0)

        # Perform PCA
        mat = Matrix(widedf[!, exps])
        mat[mat .< 1.0] .= 1.0
        mat = log.(mat)
        M = fit(PCA, mat; maxoutdim = 2)
        recon = reconstruct(M, MultivariateStats.transform(M, mat))
        error = ((recon .- mat) .^ 2)

        # Impute by SVD
        matmiss = convert(Array{Union{Float64, Missing}}, mat)
        matmiss[error .> quantile(reshape(error, :), [cutoff])] .= missing
        Impute.impute!(matmiss, Impute.SVD())
        widedf[!, exps] .= exp.(matmiss)

        ndf = stack(widedf, exps)
        rename!(ndf, "variable" => "Experiment")
        rename!(ndf, "value" => "Value")

        append!(retdf, ndf)
    end
    return sort!(retdf, ["Valency", "Cell", "subclass_1", "subclass_2", "Experiment", "%_2"])
end


function PCA_dimred()
    df = loadMixData()
    mdf = unstack(df, ["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"], "Experiment", "Value")
    mat = Matrix(mdf[!, Not(["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"])])
    Impute.impute!(mat, Impute.SVD())

    # to change Matrix type without Missing
    rmat = zeros(size(mat))
    rmat .= mat

    M = fit(PCA, rmat; maxoutdim = 1)
    mdf = mdf[!, ["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]]
    mdf."PCA" = projection(M)[:, 1]
    mdf = predictMix(mdf)
    mdf."PCA" *= mean(mdf."Predict") / mean(mdf."PCA")
    mdf[mdf."PCA" .< 1.0, "PCA"] .= 1.0
    return mdf
end
