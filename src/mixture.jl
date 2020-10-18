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


function conversionFactor(df, Vals::Vector{Int64})
    df = copy(df)
    df[!, "NewValency"] .= 0
    @assert length(Vals) == length(unique(df."Valency"))
    for (i, f) in enumerate(sort(unique(df."Valency")))
        df[(df."Valency" .== f), "NewValency"] .= Vals[i]
    end

    # Calculate predictions
    df[!, "Predict"] .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictDFRow(df[i, :])
    end

    meansh = mean(df[!, "Predict"]) ./ mean(df[!, "Value"])
    ValConvs = Dict(
        [f => mean(df[df[!, "Valency"] .== f, "Predict"]) ./ mean(df[df[!, "Valency"] .== f, "Value"]) ./ meansh for f in sort(unique(df."Valency"))],
    )
    ExpConvs = Dict(
        [d => mean(df[df[!, "Experiment"] .== d, "Predict"]) ./ mean(df[df[!, "Experiment"] .== d, "Value"]) for d in sort(unique(df."Experiment"))],
    )

    return ValConvs, ExpConvs
end

function conversionDF(df, ValConvs, ExpConvs, Vals = nothing)
    df = copy(df)
    df[!, "Adjusted"] .= df[!, "Value"]
    for i = 1:size(df)[1]
        df[i, "Adjusted"] *= ValConvs[df[i, "Valency"]]
        df[i, "Adjusted"] *= ExpConvs[df[i, "Experiment"]]
    end
    if Vals != nothing
        df[!, "NewValency"] .= 0
        @assert length(Vals) == length(unique(df."Valency"))
        for (i, f) in enumerate(sort(unique(df."Valency")))
            df[(df."Valency" .== f), "NewValency"] .= Vals[i]
        end
    end
    if !("Predict" in names(df))
        df[!, "Predict"] .= 0.0
        for i = 1:size(df)[1]
            df[i, "Predict"] = predictDFRow(df[i, :])
        end
    end
    return df
end


function fitValLoss(df, Vals)
    df = conversionDF(df, conversionFactor(df, Vals)..., Vals)
    #println(df[1, :])
    loss = 0
    for val in unique(df."Valency")
        ndf = df[df."Valency" .== val, :]
        diffs = (ndf[!, "Value"] .- ndf[!, "Predict"]) ./ ndf[!, "Predict"]
        loss += (diffs' * diffs) / length(diffs)
    end
    return loss
end


function plotValLoss(df, title = "")
    losses = zeros(10, 40)
    for val1 = 1:10
        for val2 = 1:40
            losses[val1, val2] = log(fitValLoss(df, [val1, val2]))
        end
    end
    pl = spy(losses, Guide.xlabel("Fitted valency for f=33"), Guide.ylabel("Fitted valency for f=4"), Guide.title("Log loss for " * title))
    return pl
end


function plotMixPrediction(df, title = "")
    setGadflyTheme()
    df = copy(df)
    # Calculate predictions
    df[!, "Predict"] .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictDFRow(df[i, :])
    end
    pl = plot(df, x = :Value, y = :Predict, color = :Experiment, shape = :Valency, Guide.title(title), style(key_position = :right))
    return pl
end


""" Use for only one pair of IgG subclasses and cell line / Fc receptor"""
function plotMixContinuous(df, ValConvs, ExpConvs)
    df = copy(df)
    @assert length(unique(df."Cell")) == 1
    @assert length(unique(df."subclass_1")) == 1
    @assert length(unique(df."subclass_2")) == 1
    IgGXname = Symbol(unique(df."subclass_1")[1])
    IgGYname = Symbol(unique(df."subclass_2")[1])
    # Calculate predictions
    df[!, "Predict"] .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictDFRow(df[i, :])
    end

    x = 0:0.01:1
    df4 = df[(df."Valency" .== 4), :]
    preds4 = [predictMix(df4[1, :], IgGXname, IgGYname, i, 1 - i) for i in x]
    df33 = df[(df."Valency" .== 33), :]
    preds33 = [predictMix(df33[1, :], IgGXname, IgGYname, i, 1 - i) for i in x]

    if !("Adjusted" in names(df))
        df = conversionDF(df, ValConvs, ExpConvs)
    end

    pl = plot(
        layer(x = x, y = preds4, Geom.line, Theme(default_color = colorant"red", line_width = 2px)),
        layer(x = x, y = preds33, Geom.line, Theme(default_color = colorant"green", line_width = 2px)),
        layer(df, x = "%_1", y = "Adjusted", color = "Experiment", shape = "Valency"),
        Scale.x_continuous(labels = n -> "$IgGXname $(n*100)%\n$IgGYname $(100-n*100)%"),
        Guide.xlabel(""),
        Guide.ylabel("Lbound", orientation = :vertical),
        Guide.manual_color_key("Predictions", ["f = 4", "f = 33"], ["red", "green"]),
        Guide.title("$IgGXname-$IgGYname in $(df[1, "Cell"])"),
    )

    return pl
end


function plotMixtures()
    df = loadMixData()
    draw(SVG("figure_mixture.svg", 600px, 400px), plotGrid((1, 1), [plotMixPrediction(df)]))
    draw(SVG("figure_mixture_valencies.svg", 700px, 300px), plotGrid((1, 1), [plotValLoss(df, "all data")]))
    cells = unique(df."Cell")
    pairs = unique([r."subclass_1" * "-" * r."subclass_2" for r in eachrow(df)])

    ValConvs, ExpConvs = conversionFactor(df, [4, 33])
    pls = Matrix(undef, length(cells), length(pairs))
    for (i, pair) in enumerate(pairs)
        for (j, cell) in enumerate(cells)
            pls[j, i] = plotMixPrediction(
                df[(df."Cell" .== cell) .& (df."subclass_1" .== split(pair, "-")[1]) .& (df."subclass_2" .== split(pair, "-")[2]), :],
                pair * " in " * cell,
            )
        end
    end
    draw(SVG("figure_mixture_split.svg", 2500px, 1000px), plotGrid(size(pls), pls))

    pls = Matrix(undef, length(cells), length(pairs))
    for (i, pair) in enumerate(pairs)
        for (j, cell) in enumerate(cells)
            ndf = df[(df."Cell" .== cell) .& (df."subclass_1" .== split(pair, "-")[1]) .& (df."subclass_2" .== split(pair, "-")[2]), :]
            pls[j, i] = plotMixContinuous(ndf, ValConvs, ExpConvs)
        end
    end
    draw(SVG("figure_mixture_continuous.svg", 2500px, 1000px), plotGrid(size(pls), pls))
end
