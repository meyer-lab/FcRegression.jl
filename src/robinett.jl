const RobMeasuredRecepExp = Dict(
    "FcgRI" => 232871.607,
    "FcgRIIA-131H" => 1605371.923,
    "FcgRIIA-131R" => 318818.901,
    "FcgRIIB-232I" => 394556.044,
    "FcgRIIIA-158F" => 4677645.042,
    "FcgRIIIA-158V" => 3680707.938,
)  # geometric mean

function figure2d()
    df = averageData(loadMixData("robinett/Luxetal2013-Fig2BRef.csv"))
    df = predictMix(df; recepExp = RobMeasuredRecepExp)
    df."Adjusted" = df."Value" .* (geocmean(df."Predict") / geocmean(df."Value"))
    draw(SVG("figure2d.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df)]))
end

function importRobinett()
    df = CSV.File(joinpath(FcRegression.dataDir, "robinett/Luxetal2013-Fig2Bmod.csv"), delim = ",", comment = "#") |> DataFrame
    df = dropmissing(stack(df, Not(["Cell", "Antibody", "Valency"])))
    rename!(df, ["variable" => "Experiment", "value" => "Value"])
    rename!(df, ["Antibody" => "subclass_1"])
    df[!, "%_1"] .= 1.0
    df[!, "subclass_2"] .= "None"
    df[!, "%_2"] .= 0.0
    return sort!(df, ["Valency", "Cell", "subclass_1", "subclass_2", "Experiment", "%_2"])
end

function plotRobinettCV()
    df = averageMixData(importRobinett())
    fitted = MAPLikelihood(df; robinett=true)

    # _, unfitted = fitMixMaster(df; fitKav = false, recepExp = Rtot, show_trace = false)
    # pl1 = plotPredvsMeasured(unfitted; xx = "Value", title = "With reported Kav", R2pos = (3, 1))
    pl2 = plotPredvsMeasured(fitted; xx = "Value", title = "With newly fitted Kav", R2pos = (3, 1))
    draw(SVG("figure2rob.svg", 9inch, 4inch), plotGrid((1, 2), [nothing, pl2]; sublabels = [1 1]))
end
