const RobMeasuredRecepExp = Dict(
    "FcgRI" => 232871.607, 
    "FcgRIIA-131H" => 1605371.923, 
    "FcgRIIA-131R" => 318818.901, 
    "FcgRIIB-232I" => 394556.044, 
    "FcgRIIIA-158F" => 4677645.042, 
    "FcgRIIIA-158V" => 3680707.938, 
)  # geometric mean

function figure2d()
    df = loadMixData("robinett/Luxetal2013-Fig2BRef.csv", avg = true)
    df = predictMix(df; recepExp = RobMeasuredRecepExp)
    df."Adjusted" = df."Value" .* (geocmean(df."Predict") / geocmean(df."Value"))
    draw(SVG("figure2d.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df)]))
end
