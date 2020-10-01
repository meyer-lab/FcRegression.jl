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

    return df
end


function plotMixPrediction()
    df = loadMixData()
    # Calculate predictions
    df[!, "Predict"] .= 0.0
    for i in 1:size(df)[1]
        df[i, "Predict"] = predictMixture(df[i, :])
    end

    setGadflyTheme()
    pl = plot(
        df,
        x = :Value,
        y = :Predict,
        color = :Experiment,
        shape = :Valency,
        Guide.title("Predicted vs measured mixture Lbound"),
        style(key_position = :right),
    )
    draw(SVG("figure_mixture.svg", 500px, 400px), plotGrid((1, 1), [pl]))
    return pl
end
