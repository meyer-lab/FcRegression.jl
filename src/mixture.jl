function predictMixture(dfrow)
    IgGC = zeros(size(humanIgG))
    IgGC[Symbol(dfrow."subclass_1") .== humanIgG] .= dfrow."%_1"
    IgGC[Symbol(dfrow."subclass_2") .== humanIgG] .= dfrow."%_2"

    recepExp = Dict("hFcgRIIA-131His" => 445141, "hFcgRIIB" => 31451, "hFcgRIIIA-131Val" => 657219) # geometric mean
    recepName = Dict("hFcgRIIA-131His" => "FcgRIIA-131H", "hFcgRIIB" => "FcgRIIB-232I", "hFcgRIIIA-131Val" => "FcgRIIIA-158V")
    Kav = importKav(; murine = false, retdf = true)
    Kav = Matrix(Kav[!, [recepName[dfrow."Cell"]]])
    return polyfc(1e-9, KxConst, dfrow."Valency", [recepExp[dfrow."Cell"]], IgGC, Kav).Lbound
end


function loadMixData()
    df = CSV.File(joinpath(dataDir, "lux-mixture.csv"), comment = "#") |> DataFrame!
    df = stack(df, 7:size(df)[2])
    df = dropmissing(df)
    rename!(df, "variable" => "Experiment")
    rename!(df, "value" => "Value")

    # Normalize each experiment
    for exp in unique(df.Experiment)
        meanval = mean(df[df.Experiment .== exp, "Value"])
        df[df.Experiment .== exp, "Value"] .= floor.(df[df.Experiment .== exp, "Value"] .* 1000 / meanval)
    end

    # Calculate predictions
    df[!, "Predict"] .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictMixture(df[i, :])
    end

    return df
end


function plotMixPrediction(df, title = "")
    setGadflyTheme()
    pl = plot(
        df,
        x = :Value,
        y = :Predict,
        color = :Experiment,
        shape = :Valency,
        Guide.title(title),
        style(key_position = :right),
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
            pls[j, i] = plotMixPrediction(df[(df."Cell" .== cell) .&
                                        (df."subclass_1" .== split(pair, "-")[1]) .&
                                        (df."subclass_2" .== split(pair, "-")[2]), :],
                                        pair * " in " * cell)
        end
    end
    draw(SVG("figure_mixture_split.svg", 2500px, 1000px), plotGrid(size(pls), pls))
end
