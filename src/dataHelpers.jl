using DataFrames
using Memoize
import CSV
import Distributions: fit_mle, LogNormal, pdf

const KxConst = 6.31e-13 # 10^(-12.2)

function geocmean(x)
    x = convert(Vector, x)
    x[x .<= 1.0] .= 1.0
    return geomean(x)
end

function geocstd(x)
    x = convert(Vector, x)
    x[x .<= 1.0] .= 1.0
    if length(x) <= 1
        return 0.0
    end
    return exp(std(log.(x)))
end

const murineCellTypes = ["ncMO", "cMO", "NKs", "Neu", "EO", "Kupffer", "KupfferHi"]
const humanCellTypes = ["ncMO", "cMO", "NKs", "Neu", "EO"]
const murineIgG = ["IgG1", "IgG2a", "IgG2b", "IgG3"]
const murineIgGFucose = ["IgG1", "IgG2a", "IgG2b", "IgG3", "IgG2bFucose"]
const humanIgG = ["IgG1", "IgG2", "IgG3", "IgG4"]
const murineFcgR = ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]
const humanFcgR =
    ["FcgRI", "FcgRIIA-131H", "FcgRIIA-131R", "FcgRIIB-232I", "FcgRIIB-232T", "FcgRIIC-13N", "FcgRIIIA-158F", "FcgRIIIA-158V", "FcgRIIIB"]
const humanFcgRiv = ["FcgRI", "FcgRIIA-131H", "FcgRIIA-131R", "FcgRIIB-232I", "FcgRIIIA-158F", "FcgRIIIA-158V"]
const murineActI = [1.0, -1, 1, 1]
const humanActI = [1.0, 1, 1, -1, -1, 1, 1, 1, 1]
const murineActYmax = [8e4, 5e3, 2.5e-1, 7e3, 3] # ymax for synergy plots
const humanActYmax = [5.5e4, 1.5e5, 4.5e4, 3.5e4, 3e3] # ymax for synergy plots
const dataDir = joinpath(dirname(pathof(FcRegression)), "..", "data")

function importRtot(; murine = true, genotype = "HIV", retdf = false, cellTypes = nothing)
    if murine
        df = CSV.File(joinpath(dataDir, "murine-FcgR-abundance.csv"), comment = "#") |> DataFrame
    else
        df = CSV.File(joinpath(dataDir, "human-FcgR-abundance.csv"), comment = "#") |> DataFrame
    end
    if cellTypes == nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    end
    df = combine(groupby(df, ["Cells", "Receptor"]), names(df, "Count") .=> geocmean)
    df = unstack(df, "Receptor", "Cells", "Count_geocmean")
    df = coalesce.(df, 1.0)

    if murine
        df = df[in(murineFcgR).(df.Receptor), :]
    else
        df[df[:, "Receptor"] .== "FcgRIIC", "Receptor"] .= "FcgRIIC-13N"

        generic_type = ["FcgRIIA", "FcgRIIB", "FcgRIIIA"]
        prefixes = ["FcgRIIA-131", "FcgRIIB-232", "FcgRIIIA-158"]
        options = [['H', 'R'], ['I', 'T'], ['F', 'V']]
        ncols = size(df)[2] - 1

        for i = 1:3
            rowidx = findfirst(df[:, "Receptor"] .== generic_type[i])
            if genotype[i] == options[i][1]
                df[rowidx, "Receptor"] = prefixes[i] * options[i][1]
                push!(df, [prefixes[i] * options[i][2]; repeat([0.0], ncols)])
            elseif genotype[i] == options[i][2]
                df[rowidx, "Receptor"] = prefixes[i] * options[i][2]
                push!(df, [prefixes[i] * options[i][1]; repeat([0.0], ncols)])
            else  # heterozygous
                push!(df, [prefixes[i] * options[i][1]; Array(df[rowidx, 2:end]) ./ 2])
                push!(df, [prefixes[i] * options[i][2]; Array(df[rowidx, 2:end]) ./ 2])
                df = df[df[:, "Receptor"] .!= generic_type[i], :]
            end
        end
        sort!(df, ["Receptor"])
    end
    @assert df.Receptor == (murine ? murineFcgR : humanFcgR)
    if retdf
        return deepcopy(df[!, ["Receptor"; names(df)[in(cellTypes).(names(df))]]])
    else
        return deepcopy(Matrix{Float64}(df[!, cellTypes]))
    end
end


""" Import human or murine affinity data. """
@memoize function importKav(; murine = true, c1q = false, invitro = false, IgG2bFucose = false, retdf = false)
    if murine
        df = CSV.File(joinpath(dataDir, "murine-affinities.csv"), comment = "#") |> DataFrame
    else
        df = CSV.File(joinpath(dataDir, "human-affinities.csv"), comment = "#") |> DataFrame
    end

    IgGlist = copy(murine ? murineIgG : humanIgG)
    FcRecep = copy(murine ? murineFcgR : (invitro ? humanFcgRiv : humanFcgR))
    if IgG2bFucose
        append!(IgGlist, ["IgG2bFucose"])
    end
    if c1q
        append!(FcRecep, ["C1q"])
    end
    df = stack(df; variable_name = "IgG", value_name = "Kav")
    df = unstack(df, "FcgR", "Kav")
    dropmissing!(df)
    df = df[in(IgGlist).(df.IgG), :]

    if retdf
        return df[!, ["IgG"; FcRecep]]
    else
        return Matrix{Float64}(df[!, FcRecep])
    end
end


""" Import cell depletion data. """
function importDepletion(dataType)
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

    affinity = importKav(murine = true, c1q = c1q, IgG2bFucose = true, retdf = true)
    df = leftjoin(df, affinity, on = "Condition" => "IgG")

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


""" Import measurements of receptor amounts. """
@memoize function importInVitroRtotDist(robinett = false)
    local df
    if robinett
        df = CSV.File(joinpath(FcRegression.dataDir, "robinett/FcgRquant.csv"), delim = ",", comment = "#") |> DataFrame
    else
        df = CSV.File(joinpath(dataDir, "receptor_amount_mar2021.csv"), delim = ",", comment = "#") |> DataFrame
    end
    return [fit_mle(LogNormal, (df[df."Receptor" .== rcp, "Measurements"])) for rcp in humanFcgRiv]
end


@memoize function importKavDist(; inflation = 0.1, retdf = true)
    df = CSV.File(joinpath(dataDir, "FcgR-Ka-Bruhns_with_variance.csv"), delim = ",", comment = "#") |> DataFrame
    function parstr(x)
        params = parse.(Float64, split(x, "|"))
        params .+= inflation
        params .*= 1e5      # Bruhns data is written in 1e5 units
        return LogNormal(log(params[1]), 2.3) # std of 10x
    end
    xdf = parstr.(df[:, Not("IgG")])
    insertcols!(xdf, 1, "IgG" => df[:, "IgG"])
    if retdf
        return xdf
    else
        return Matrix(xdf[:, Not("IgG")])
    end
end

const f4Dist = LogNormal(log(4), 0.1)
const f33Dist = LogNormal(log(33), 0.1)
const KxStarDist = LogNormal(log(KxConst), 2.0)   # ~ 4.37 in Robinett
