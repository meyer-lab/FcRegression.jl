""" Import cell depletion data. """
function importDepletion(dataType; Kav::Union{AbstractDataFrame, Nothing} = nothing)
    c1q = false
    if dataType == "ITP"
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
    df[!, "Target"] = 1.0 .- df[!, "Target"] ./ 100.0
    if "Neutralization" in names(df)
        neut = -log.(df[!, "Neutralization"] / 50.0)
        df[!, "Neutralization"] .= replace!(neut, Inf => 0.0)
    end

    Kavd = importKav(murine = true, c1q = c1q, IgG2bFucose = true, retdf = true)
    if Kav !== nothing
        Kav[Kav."IgG" .== "IgG2c", "IgG"] .= "IgG2a"
        # replace tha value in
        for igg in Kav."IgG"
            Kavd[Kavd."IgG" .== igg, names(Kav)[2:end]] = Kav[Kav."IgG" .== igg, names(Kav)[2:end]]
        end
    end
    df = leftjoin(df, Kavd, on = "Condition" => "IgG")

    # The mG053 antibody doesn't bind to the virus
    if dataType == "HIV"
        df[df[:, "Label"] .== "mG053", ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]] .= 0.0
    end

    df[df[:, "Background"] .== "R1KO", "FcgRI"] .= 0.0
    df[df[:, "Background"] .== "R2KO", "FcgRIIB"] .= 0.0
    df[df[:, "Background"] .== "R3KO", "FcgRIII"] .= 0.0
    df[df[:, "Background"] .== "R1/3KO", ["FcgRI", "FcgRIII"]] .= 0.0
    df[df[:, "Background"] .== "R1/4KO", ["FcgRI", "FcgRIV"]] .= 0.0
    df[df[:, "Background"] .== "R4block", "FcgRIV"] .= 0.0
    df[df[:, "Background"] .== "gcKO", ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]] .= 0.0
    df[df[:, "Condition"] .== "IgG1D265A", ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]] .= 0.0

    for pair in ["R" => "FcγR", "1" => "I", "2" => "II", "3" => "III", "4" => "IV", "gc" => "γc"]
        df[!, "Background"] = map(x -> replace(x, pair), df.Background)
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
    elseif dataType == "ITP"
        df = CSV.File(joinpath(dataDir, "schwab_ITP_humanized.csv"), delim = ",", comment = "#") |> DataFrame
        df = stack(df, ["IgG1", "IgG2", "IgG3", "IgG4"])
        df = disallowmissing!(df[completecases(df), :])
        rename!(df, ["variable" => "Condition", "value" => "Target"])

        df[!, "Target"] .= 1.0 .- df.Target ./ 100.0
        affinity = importKav(murine = false, c1q = false, retdf = true)
    else
        @error "Data type not found"
    end

    df = leftjoin(df, affinity, on = "Condition" => "IgG")
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
    return df
end
