function predictMix(dfrow, IgGXname, IgGYname, IgGX, IgGY)
    IgGC = zeros(size(humanIgG))
    IgGC[IgGXname .== humanIgG] .= IgGX
    IgGC[IgGYname .== humanIgG] .= IgGY

    recepExp = Dict("hFcgRIIA-131His" => 445141, "hFcgRIIB" => 31451, "hFcgRIIIA-131Val" => 657219)  # geometric mean
    recepName = Dict("hFcgRIIA-131His" => "FcgRIIA-131H", "hFcgRIIB" => "FcgRIIB-232I", "hFcgRIIIA-131Val" => "FcgRIIIA-158V")
    Kav = importKav(; murine = false, retdf = true)
    Kav = Matrix(Kav[!, [recepName[dfrow."Cell"]]])
    return polyfc(1e-9, KxConst, dfrow."Valency", [recepExp[dfrow."Cell"]], IgGC, Kav).Lbound
end


function predictDFRow(dfrow)
    return predictMix(dfrow, Symbol(dfrow."subclass_1"), Symbol(dfrow."subclass_2"), dfrow."%_1", dfrow."%_2")
end


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


function plotMixPrediction(df, title = "")
    setGadflyTheme()
    # Calculate predictions
    df[!, "Predict"] .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictDFRow(df[i, :])
    end
    pl = plot(df, x = :Value, y = :Predict, color = :Experiment, shape = :Valency, Guide.title(title), style(key_position = :right))
    return pl
end


""" Use for only one pair of IgG subclasses and cell line / Fc receptor"""
function plotMixContinuous(df; split_val = true)
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
    if split_val
        for val in unique(df."Valency")
            for exp in unique(df."Experiment")
                meanval = mean(df[(df."Experiment" .== exp) .& (df."Valency" .== val), "Value"])
                meanpred = mean(df[(df."Experiment" .== exp) .& (df."Valency" .== val), "Predict"])
                df[(df."Experiment" .== exp) .& (df."Valency" .== val), "Value"] .*= meanpred / meanval
            end
        end
    else
        for exp in unique(df."Experiment")
            meanval = mean(df[df."Experiment" .== exp, "Value"])
            meanpred = mean(df[df."Experiment" .== exp, "Predict"])
            df[df."Experiment" .== exp, "Value"] .*= meanpred / meanval
        end
    end

    pl = plot(
        layer(x = x, y = preds4, Geom.line, Theme(default_color = colorant"red", line_width = 2px)),
        layer(x = x, y = preds33, Geom.line, Theme(default_color = colorant"green", line_width = 2px)),
        layer(df, x = "%_1", y = "Value", color = "Experiment", shape = "Valency"),
        Scale.x_continuous(labels = n -> "$IgGXname $(n*100)%\n$IgGYname $(100-n*100)%"),
        Guide.xlabel(""),
        Guide.ylabel("Lbound", orientation = :vertical),
        Guide.manual_color_key("Predictions", ["f = 4", "f = 33"], ["red", "green"]),
        Guide.title("$IgGXname-$IgGYname in $(df[1, "Cell"])"),
    )

    return pl
end


function fitValLoss(df, val1, val2)
    df = copy(df)
    @assert length(unique(df."Cell")) == 1
    @assert length(unique(df."subclass_1")) == 1
    @assert length(unique(df."subclass_2")) == 1
    IgGXname = Symbol(unique(df."subclass_1")[1])
    IgGYname = Symbol(unique(df."subclass_2")[1])

    rename!(df, "Valency" => "OldValency")
    df[!, "Valency"] .= 1
    df[(df."OldValency" .== 4), "Valency"] .= val1
    df[(df."OldValency" .== 33), "Valency"] .= val2
    # Calculate predictions
    df[!, "Predict"] .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictDFRow(df[i, :])
    end

    for exp in unique(df."Experiment")
        meanval = mean(df[(df."Experiment" .== exp), "Value"])
        meanpred = mean(df[(df."Experiment" .== exp), "Predict"])
        df[(df."Experiment" .== exp), "Value"] .*= meanpred / meanval
    end

    loss = 0
    for val in unique(df."Valency")
        ndf = df[df."Valency" .== val, :]
        diffs = (ndf[!, "Value"] .- ndf[!, "Predict"]) ./ ndf[!, "Predict"]
        loss += (diffs' * diffs) / length(diffs)
    end
    return loss
end


function plotValLoss(df, title = "")
    losses = zeros(10, 33)
    for val1 in 1:10
        for val2 in 1:33
            losses[val1, val2] = log(fitValLoss(df, val1, val2))
        end
    end
    pl = spy(
        losses,
        Guide.xlabel("Fitted valency for f=33"),
        Guide.ylabel("Fitted valency for f=4"),
        Guide.title("Log loss for " * title),
    )
    return pl
end



function plotMixtures()
    df = loadMixData()
    draw(SVG("figure_mixture.svg", 600px, 400px), plotGrid((1, 1), [plotMixPrediction(df)]))
    cells = unique(df."Cell")
    pairs = unique([r."subclass_1" * "-" * r."subclass_2" for r in eachrow(df)])
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

    plson = Matrix(undef, length(cells), length(pairs))
    plsoff = Matrix(undef, length(cells), length(pairs))
    plv = Matrix(undef, length(cells), length(pairs))
    for (i, pair) in enumerate(pairs)
        for (j, cell) in enumerate(cells)
            ndf = df[(df."Cell" .== cell) .& (df."subclass_1" .== split(pair, "-")[1]) .& (df."subclass_2" .== split(pair, "-")[2]), :]
            plson[j, i] = plotMixContinuous(ndf)
            plsoff[j, i] = plotMixContinuous(ndf; split_val = false)
            plv[j, i] = plotValLoss(ndf, pair * " in " * cell)
        end
    end
    draw(SVG("figure_mixture_continuous.svg", 2500px, 1000px), plotGrid(size(plson), plson))
    draw(SVG("figure_mixture_continuous_notsplit.svg", 2500px, 1000px), plotGrid(size(plsoff), plsoff))
    draw(SVG("figure_mixture_valencies.svg", 2500px, 1000px), plotGrid(size(plv), plv))
end
