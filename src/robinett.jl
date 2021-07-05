const RobMeasuredRecepExp = Dict(
    "FcgRI" => 232871.607,
    "FcgRIIA-131H" => 1605371.923,
    "FcgRIIA-131R" => 318818.901,
    "FcgRIIB-232I" => 394556.044,
    "FcgRIIIA-158F" => 4677645.042,
    "FcgRIIIA-158V" => 3680707.938,
)  # geometric mean

RobFitMeasuredRecepExp_avg = Dict(
    "FcgRI" => 11593.99461864337,
    "FcgRIIA-131H" => 79926.76168625052,
    "FcgRIIA-131R" => 61350.26603001091,
    "FcgRIIB-232I" => 277772.301710943,
    "FcgRIIIA-158F" => 232886.23350665424,
    "FcgRIIIA-158V" => 183251.65775134592,
)
RobFitMeasuredRecepExp = Dict(
    "FcgRI" => 11593.99461864337,
    "FcgRIIA-131H" => 93102.55658789474,
    "FcgRIIA-131R" => 78843.09186661601,
    "FcgRIIB-232I" => 356908.7172811272,
    "FcgRIIIA-158F" => 232886.23350665424,
    "FcgRIIIA-158V" => 214391.89886947838,
)

RobFit_valencies_avg = [1.8412288223384095, 12.140021558657603]
RobFit_valencies = [2.3132889097502374, 12.140021558657603]

RobFitKx_avg = 5.753975202649211e-16
RobFitKx = 5.753975202649498e-16

RobConversion_facs = [1.7411125324367682, 1.1851999547787864, 0.7366318246623974, 2.8302946874787773]

function figure2d(avg = false)
    df = avg ? averageData(loadMixData("robinett/Luxetal2013-Fig2BRef.csv")) : loadMixData("robinett/Luxetal2013-Fig2BRef.csv")

    df."NewValency" = convert.(Float64, df."Valency")
    df[df[!, "NewValency"] .== 4.0, "NewValency"] .= avg ? RobFit_valencies_avg[1] : RobFit_valencies[1]
    df[(df[!, "NewValency"]) .== 26.0, "NewValency"] .= avg ? RobFit_valencies_avg[2] : RobFit_valencies[2]

    df = predictMix(df; recepExp = avg ? RobFitMeasuredRecepExp_avg : RobFitMeasuredRecepExp, KxStar = avg ? RobFitKx_avg : RobFitKx)
    df[!, :Adjusted] = df[!, "Value"]
    if !avg
        df = fit_experiments(RobConversion_facs, df)
    end
    #df."Adjusted" = df."Value" .* (geocmean(df."Predict") / geocmean(df."Value"))
    draw(SVG("figure2d.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df)]))
end
