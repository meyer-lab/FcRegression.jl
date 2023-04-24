""" Humanized mice data from Lux 2014, Schwab 2015 """
function importHumanized(dataType)
    if dataType == "ITP"
        df = CSV.File(joinpath(dataDir, "schwab_ITP_humanized.csv"), delim = ",", comment = "#") |> DataFrame
        df = stack(df, ["IgG1", "IgG2", "IgG3", "IgG4"])
        df = disallowmissing!(df[completecases(df), :])
        rename!(df, ["variable" => "Condition", "value" => "Target"])
        df[!, "Target"] .= 1.0 .- df.Target ./ 100.0
        df."Genotype" = [g[1] * "I" * g[3] for g in df."Genotype"]   # FcgRIIB default as 232I
        @error "Data type not found"
    end
    @assert maximum(df."Target") <= 1.0
    return df
end
