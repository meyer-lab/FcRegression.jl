using MultivariateStats
using StatsBase

""" Load mixture in vitro binding data """
function loadMixData(fn = "CHO_IgG_mixture_binding.csv")
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame

    df = stack(df, Not(["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])
    df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0

    baseline = combine(groupby(df, "Experiment"), "Value" => geomean => "Baseline")
    df = innerjoin(df, baseline, on = "Experiment")
    df[!, "Value"] ./= df[!, "Baseline"]    # normalize fluorescence by daily geomean
    df = df[!, Not(["Experiment", "Baseline"])]
    rename!(df, "Cell" => "Receptor")

    df[!, "%_1"] ./= 100.0
    df[!, "%_2"] ./= 100.0

    replace!(df."Receptor", "CHO-hFcgRIIA-131His" => "FcgRIIA-131H")
    replace!(df."Receptor", "CHO-hFcgRIIB" => "FcgRIIB-232I")
    replace!(df."Receptor", "CHO-hFcgRIIIA-131Val" => "FcgRIIIA-158V")
    replace!(df."Receptor", "CHO-FcgRIA" => "FcgRI")
    replace!(df."Receptor", "CHO-hFcgRIIA-131Arg" => "FcgRIIA-131R")
    replace!(df."Receptor", "CHO-hFcgRIIIA-158Phe" => "FcgRIIIA-158F")

    return sort!(df, ["Valency", "Receptor", "subclass_1", "subclass_2", "%_2"])
end

function importRobinett()
    df = CSV.File(joinpath(dataDir, "robinett_binding.csv"), delim = ",", comment = "#") |> DataFrame
    for i = 1:4
        cn = "Replicate $i"
        df[!, cn] ./= geomean(df[Not(ismissing.(df[!, cn])), cn])
    end
    df = dropmissing(stack(df, Not(["Receptor", "Antibody", "Valency"])))
    rename!(df, ["variable" => "Experiment", "value" => "Value"])
    rename!(df, ["Antibody" => "Subclass"])

    return sort!(df, ["Valency", "Receptor", "Subclass", "Experiment"])
end


""" Make statistics of individual cell types and subclass types """
function averageMixData(df = loadMixData(); combSingle = false)
    # Combine cases of single IgGs into one entry
    if combSingle
        df[df."%_1" .== 0.0, "subclass_1"] .= "None"
        df[df."subclass_1" .== "None", "%_1"] .= 1.0
        df[df."subclass_1" .== "None", "%_2"] .= 0.0
        df[df."subclass_1" .== "None", "subclass_1"] = df[df."subclass_1" .== "None", "subclass_2"]
        df[df."%_2" .== 0.0, "subclass_2"] .= "None"
    end

    valname = "Adjusted" in names(df) ? "Adjusted" : "Value"
    combining = ["Value"]
    if "Experiment" in names(df)
        append!(combining, ["Experiment"])
    end
    return combine(groupby(df, Not(combining)), valname => StatsBase.median => valname, valname => lower => "xmin", valname => upper => "xmax")
end

""" Transform combined single IgG case back to pairs"""
function combSing2pair(df)
    @assert !("None" in df."subclass_1")
    @assert "None" in df."subclass_2"

    ndf = copy(df[[], :])
    isotypes = unique(df."subclass_1")
    for row in eachrow(df)
        if row."subclass_2" == "None"
            for igg in isotypes
                crow = deepcopy(row)
                if igg > row."subclass_1"
                    crow."subclass_2" = igg
                    push!(ndf, crow)
                elseif igg < row."subclass_1"
                    crow."subclass_2" = row."subclass_1"
                    crow."subclass_1" = igg
                    crow."%_2" = 1.0
                    crow."%_1" = 0.0
                    push!(ndf, crow)
                end
            end
        else
            push!(ndf, row)
        end
    end
    return sort!(ndf, names(df)[in(["Valency", "Receptor", "subclass_1", "subclass_2", "Experiment", "%_2"]).(names(df))])
end
