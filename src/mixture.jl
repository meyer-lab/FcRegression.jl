function loadMixData()
    df = CSV.File(joinpath(dataDir, "lux-mixture.csv"), comment = "#") |> DataFrame!
    df = stack(df, 7:size(df)[2])
    df = dropmissing(df)
    rename!(df, "variable" => "Experiment")
    rename!(df, "value" => "Value")
    df[!, "Value"] = convert.(Float64, df[!, "Value"])

    df[!, "%_1"] ./= 100.0
    df[!, "%_2"] ./= 100.0

    return df
end


function predictMix(dfrow, IgGXname, IgGYname, IgGX, IgGY)
    IgGC = zeros(size(humanIgG))
    IgGC[IgGXname .== humanIgG] .= IgGX
    IgGC[IgGYname .== humanIgG] .= IgGY

    recepExp = Dict("hFcgRIIA-131His" => 445141, "hFcgRIIB" => 31451, "hFcgRIIIA-131Val" => 657219)  # geometric mean
    recepName = Dict("hFcgRIIA-131His" => "FcgRIIA-131H", "hFcgRIIB" => "FcgRIIB-232I", "hFcgRIIIA-131Val" => "FcgRIIIA-158V")
    Kav = importKav(; murine = false, retdf = true)
    Kav = Matrix(Kav[!, [recepName[dfrow."Cell"]]])
    val = "NewValency" in names(dfrow) ? dfrow."NewValency" : dfrow."Valency"
    return polyfc(1e-9, KxConst, val, [recepExp[dfrow."Cell"]], IgGC, Kav).Lbound
end


function predictDFRow(dfrow)
    return predictMix(dfrow, Symbol(dfrow."subclass_1"), Symbol(dfrow."subclass_2"), dfrow."%_1", dfrow."%_2")
end


function TwoDFit(X::Matrix, Y::Matrix)
    """
    Fit X_ij * p_i * q_j ≈ Y_ij
    p[1] == 1 (fixed)
    length(p) == size(X, 1) - 1 == m
    length(q) == size(X, 2) == n
    v == vcat(p, q)
    """
    @assert size(X) == size(Y)
    m, n = size(X)
    f(p::Vector, q::Vector) = sum((reshape([1.0; p], :, 1) .* X .* reshape(q, 1, :) .- Y) .^ 2)
    f(v::Vector) = f(v[1:(m - 1)], v[m:end])
    init_v = ones(m + n - 1)
    od = OnceDifferentiable(f, init_v; autodiff = :forward)
    res = optimize(od, init_v, BFGS()).minimizer
    return [1.0; res[1:(m - 1)]], res[m:end]
end


function MixtureFitLoss(df, ValConv::Vector, ExpConv::Vector; logscale = false)
    df = copy(df)
    df[!, "Predict"] .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictDFRow(df[i, :])
    end
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
    nv, ne = length(unique(df."Valency")), length(unique(df."Experiment"))
    f(p::Vector, q::Vector) = MixtureFitLoss(df, [1.0; p], q; logscale = logscale)[1]
    f(v::Vector) = f(v[1:(nv - 1)], v[nv:end])
    init_v = ones(nv + ne - 1)
    od = OnceDifferentiable(f, init_v; autodiff = :forward)
    res = optimize(od, init_v, BFGS()).minimizer
    p, q = [1.0; res[1:(nv - 1)]], res[nv:end]
    res = MixtureFitLoss(df, p, q; logscale = logscale)
    return Dict("loss" => res[1], "df" => res[2], "ValConv" => p, "ExpConv" => q)
end


function plotValLoss(df, title = ""; logscale = false)
    df = copy(df)
    df[!, "NewValency"] .= 0

    losses = zeros(10, 40)
    for val1 = 1:10
        df[(df."Valency" .== 4), "NewValency"] .= val1
        for val2 = 1:40
            df[(df."Valency" .== 33), "NewValency"] .= val2
            losses[val1, val2] = log(MixtureFit(df; logscale = logscale)["loss"])
        end
    end
    pl = spy(losses, Guide.xlabel("Fitted valency for f=33"), Guide.ylabel("Fitted valency for f=4"), Guide.title("Log loss for " * title))
    return pl
end


function plotMixPrediction(df, title = "")
    @assert "Predict" in names(df)
    pl = plot(df, x = :Value, y = :Predict, color = :Experiment, shape = :Valency, Guide.title(title), style(key_position = :right))
    return pl
end


""" Use for only one pair of IgG subclasses and cell line / Fc receptor"""
function plotMixContinuous(df; logscale = false)
    df = copy(df)
    @assert length(unique(df."Cell")) == 1
    @assert length(unique(df."subclass_1")) == 1
    @assert length(unique(df."subclass_2")) == 1
    IgGXname = Symbol(unique(df."subclass_1")[1])
    IgGYname = Symbol(unique(df."subclass_2")[1])

    x = 0:0.01:1
    df4 = df[(df."Valency" .== 4), :]
    preds4 = [predictMix(df4[1, :], IgGXname, IgGYname, i, 1 - i) for i in x]
    df33 = df[(df."Valency" .== 33), :]
    preds33 = [predictMix(df33[1, :], IgGXname, IgGYname, i, 1 - i) for i in x]

    @assert "Adjusted" in names(df)

    pl = plot(
        layer(x = x, y = preds4, Geom.line, Theme(default_color = colorant"red", line_width = 2px)),
        layer(x = x, y = preds33, Geom.line, Theme(default_color = colorant"green", line_width = 2px)),
        layer(df, x = "%_1", y = "Adjusted", color = "Experiment", shape = "Valency"),
        Scale.x_continuous(labels = n -> "$IgGXname $(n*100)%\n$IgGYname $(100-n*100)%"),
        (logscale ? Scale.y_log10 : Scale.y_continuous),
        Guide.xlabel(""),
        Guide.ylabel("Lbound", orientation = :vertical),
        Guide.manual_color_key("Predictions", ["f = 4", "f = 33"], ["red", "green"]),
        Guide.title("$IgGXname-$IgGYname in $(df[1, "Cell"])"),
    )

    return pl
end


function plotMixtures()
    setGadflyTheme()
    res1 = MixtureFit(loadMixData(); logscale = false)
    println("Linear Scale: Fitted ValConv = ", res1["ValConv"], ", ExpConv = ", res1["ExpConv"])
    df1 = res1["df"]
    res2 = MixtureFit(loadMixData(); logscale = true)
    println("Log Scale: Fitted ValConv = ", res2["ValConv"], ", ExpConv = ", res2["ExpConv"])
    df2 = res2["df"]
    #draw(SVG("figure_mixture_valencies_linear.svg", 700px, 300px), plotGrid((1, 1), [plotValLoss(df1, "all data"; logscale = false)]))
    #draw(SVG("figure_mixture_valencies_log.svg", 700px, 300px), plotGrid((1, 1), [plotValLoss(df2, "all data"; logscale = true)]))
    cells = unique(df1."Cell")
    pairs = unique([r."subclass_1" * "-" * r."subclass_2" for r in eachrow(df1)])

    pls1 = Matrix(undef, length(cells), length(pairs))
    pls2 = Matrix(undef, length(cells), length(pairs))
    for (i, pair) in enumerate(pairs)
        for (j, cell) in enumerate(cells)
            ndf1 = df1[(df1."Cell" .== cell) .& (df1."subclass_1" .== split(pair, "-")[1]) .& (df1."subclass_2" .== split(pair, "-")[2]), :]
            pls1[j, i] = plotMixContinuous(ndf1; logscale = false)
            ndf2 = df2[(df2."Cell" .== cell) .& (df2."subclass_1" .== split(pair, "-")[1]) .& (df2."subclass_2" .== split(pair, "-")[2]), :]
            pls2[j, i] = plotMixContinuous(ndf2; logscale = true)
        end
    end
    draw(SVG("figure_mixture_continuous_linear.svg", 2500px, 1000px), plotGrid(size(pls1), pls1))
    draw(SVG("figure_mixture_continuous_log.svg", 2500px, 1000px), plotGrid(size(pls2), pls2))
end
