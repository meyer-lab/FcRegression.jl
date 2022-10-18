""" Import cell depletion data. """
function importDepletion(dataType)
    c1q = false
    if dataType == "melanomaNew"
        filename = "science2005-melanoma.csv"
    elseif dataType == "ITP"
        filename = "nimmerjahn-ITP.csv"
    elseif dataType == "blood"
        filename = "nimmerjahn-CD20-blood.csv"
        c1q = true
    elseif dataType == "bone"
        filename = "nimmerjahn-CD20-bone.csv"
        c1q = true
    elseif dataType == "melanoma"
        filename = "nimmerjahn-melanoma.csv"
    elseif dataType == "HIV"
        filename = "elsevier-HIV.csv"
    elseif dataType == "Bcell"
        filename = "Lux_et_al_C57BL6.csv"
        c1q = true
    else
        @error "Data type not found"
    end

    df = CSV.File(joinpath(dataDir, filename), delim = ",", comment = "#") |> DataFrame

    if "Source" in names(df)
        df = df[!, Not("Source")]
    end
    if !("Target" in names(df))  # Target: the larger, the greater effect
        df."Target" = 1 .- df."Measurement" ./ df."Baseline"
        df[df."Target" .> 1.0, "Target"] .= 1.0
        df[df."Target" .< 0.0, "Target"] .= 0.0
    else
        df[!, "Target"] = 1.0 .- df[!, "Target"] ./ 100.0
    end
    if "Subclass" in names(df)
        rename!(df, "Subclass" => "Condition")
    end
    if "Dose" in names(df)
        rename!(df, "Dose" => "Concentration")
    end
    if "Neutralization" in names(df)
        neut = -log.(df[!, "Neutralization"] / 50.0)
        df[!, "Neutralization"] .= replace!(neut, Inf => 0.0)
    end
    return df
end


""" Humanized mice data from Lux 2014, Schwab 2015 """
function importHumanized(dataType)
    if dataType in ["blood", "spleen", "bone"]
        df = CSV.File(joinpath(dataDir, "lux_humanized_CD19.csv"), delim = ",", comment = "#") |> DataFrame
        df = dropmissing(df, Symbol(dataType), disallowmissing = true)
        df[!, "Target"] = 1.0 .- df[!, Symbol(dataType)] ./ 100.0
        df[!, "Condition"] .= "IgG1"
        df = df[!, ["Genotype", "Concentration", "Condition", "Target"]]
        affinity = importKav(murine = false, c1q = true, retdf = true)
        df = leftjoin(df, affinity, on = "Condition" => "IgG")
    elseif dataType == "ITP"
        df = CSV.File(joinpath(dataDir, "schwab_ITP_humanized.csv"), delim = ",", comment = "#") |> DataFrame
        df = stack(df, ["IgG1", "IgG2", "IgG3", "IgG4"])
        df = disallowmissing!(df[completecases(df), :])
        rename!(df, ["variable" => "Condition", "value" => "Target"])
        df[!, "Target"] .= 1.0 .- df.Target ./ 100.0
        df."Genotype" = [g[1]*"I"*g[3] for g in df."Genotype"]   # FcgRIIB default as 232I
    elseif dataType == "bloodFig5c"
        df = CSV.File(joinpath(dataDir, "cell_report_2014_fig5c.csv"), delim = ",", comment = "#") |> DataFrame
        df[!, "Target"] = df."Target" ./ 100.0
        df[!, "Target"] = 1 .- df."Target"
    else
        @error "Data type not found"
    end
    @assert maximum(df."Target") <= 1.0
    return df
end


""" Mouse mix IgG depletion data from Lux """
function importDeplExp()
    df = CSV.File(joinpath(dataDir, "lux_depletion_mixedIgG_sep2021.csv"), delim = ",", comment = "#") |> DataFrame
    @assert all([item in names(df) for item in ["subclass_1", "%_1", "subclass_2", "%_2"]])
    df[ismissing.(df."%_1"), "%_1"] .= 0
    df[ismissing.(df."%_2"), "%_2"] .= 0
    df[ismissing.(df."subclass_1"), "subclass_1"] .= "PBS"
    df[ismissing.(df."subclass_2"), "subclass_2"] .= "PBS"
    df = coalesce.(df, 0)
    df."depletion" /= 100.0
    df[df."depletion" .< 0.0, "depletion"] .= 0.0
    rename!(df, "depletion" => "Target")

    df = combine(
        groupby(df, Not("Target")),
        "Target" => geomean => "Target",
        "Target" => (x -> quantile(x, 0.25)) => "xmin",
        "Target" => (x -> quantile(x, 0.75)) => "xmax",
    )
    return df
end
