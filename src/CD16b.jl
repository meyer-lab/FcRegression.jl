function importCD16b()
    df = CSV.File(joinpath(FcRegression.dataDir, "dec2022-FcgR3b-binding.csv"), delim = ",", comment = "#") |> DataFrame
    df = stack(df, Not(["Valency", "Cell", "Subclass"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])
    df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0

    expression = DataFrame(Dict("Experiment" => unique(df."Experiment"),
                            "CHO-FcgRIIIB-NA1" => [77.1, 90.8, 94.1, 90.8, 95.8, 98.0], 
                            "CHO-FcgRIIIB-NA2" => [53.3, 81.7, 91.8, 86.6, 92.4, 92.9]))
    # Expression data for "07.10.19 AM" was not available; used geomean of others to impute
    expression = stack(expression, Not("Experiment"), variable_name = "Cell", value_name = "Expression")
    df = innerjoin(df, expression, on = ["Cell", "Experiment"])

    rename!(df, "Cell" => "Receptor")
    replace!(df."Receptor", "CHO-FcgRIIIB-NA1" => "FcgRIIIB-NA1")
    replace!(df."Receptor", "CHO-FcgRIIIB-NA2" => "FcgRIIIB-NA2")
    return df
end

function inferCD16b()
    CD16bRtot = Dict("FcgRIIIB-NA1" => 100_000.0, "FcgRIIIB-NA2" => 100_000.0)
    pred = FcRegression.predMix(df; Kav = Kav[!, [1,8,9]], Rtot = CD16bRtot)
end
