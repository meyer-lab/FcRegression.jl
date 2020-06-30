using DataFrames
import CSV
import StatsBase.geomean
using Memoize

const KxConst = 6.31e-13 # 10^(-12.2)

function geocmean(x)
    x = convert(Vector, x)
    x[x .<= 1.0] .= 1.0
    return geomean(x)
end

const cellTypes = [:ncMO, :cMO, :NKs, :Neu, :EO]
const murineIgG = [:IgG1, :IgG2a, :IgG2b, :IgG3]
const humanIgG = [:IgG1, :IgG2, :IgG3, :IgG4]
const murineFcgR = [:FcgRI, :FcgRIIB, :FcgRIII, :FcgRIV]
const humanFcgR =
    Symbol.(["FcgRI", "FcgRIIA-131H", "FcgRIIA-131R", "FcgRIIB-232I", "FcgRIIB-232T", "FcgRIIC-13N", "FcgRIIIA-158V", "FcgRIIIA-158F", "FcgRIIIB"])
const murineActI = [1, -1, 1, 1]
const humanActI = [1, 1, 1, -1, -1, 1, 1, 1, 1]
const dataDir = joinpath(dirname(pathof(FcgR)), "..", "data")

@memoize function importRtot(; murine = true, genotype = "HIV", retdf = false)
    if murine
        df = CSV.File(joinpath(dataDir, "murine-FcgR-abundance.csv"), comment = "#") |> DataFrame!
    else
        df = CSV.File(joinpath(dataDir, "human-FcgR-abundance.csv"), comment = "#") |> DataFrame!
    end
    df = combine(groupby(df, [:Cells, :Receptor]), names(df, :Count) .=> geocmean)
    df = unstack(df, :Receptor, :Cells, :Count_geocmean)
    dropmissing!(df)
    df[!, :Receptor] = map(Symbol, df[!, :Receptor])
    if murine
        df = df[in(murineFcgR).(df.Receptor), :]
    else
        df[df[:, :Receptor] .== :FcgRIIC, :Receptor] .= Symbol("FcgRIIC-13N")

        generic_type = Symbol.(["FcgRIIA", "FcgRIIB", "FcgRIIIA"])
        prefixes = ["FcgRIIA-131", "FcgRIIB-232", "FcgRIIIA-158"]
        options = [['H', 'R'], ['I', 'T'], ['V', 'F']]
        ncols = size(df)[2] - 1
        for i = 1:3
            rowidx = findfirst(df[:, :Receptor] .== generic_type[i])
            if genotype[i] == options[i][1]
                df[rowidx, :Receptor] = Symbol(prefixes[i] * options[i][1])
                insert!.(df, rowidx + 1, [Symbol(prefixes[i] * options[i][2]); repeat([0.0], ncols)])
            elseif genotype[i] == options[i][2]
                df[rowidx, :Receptor] = Symbol(prefixes[i] * options[i][2])
                insert!.(df, rowidx, [Symbol(prefixes[i] * options[i][1]); repeat([0.0], ncols)])
            else  # heterozygous
                insert!.(df, rowidx, [Symbol(prefixes[i] * options[i][1]); Array(df[rowidx, 2:end]) ./ 2])
                insert!.(df, rowidx + 1, [Symbol(prefixes[i] * options[i][2]); Array(df[rowidx, 2:end])])
                df = df[df[:, :Receptor] .!= generic_type[i], :]
            end
        end

    end
    @assert df.Receptor == (murine ? murineFcgR : humanFcgR)
    if retdf
        return df[!, [:Receptor; cellTypes]]
    else
        return convert(Matrix{Float64}, df[!, cellTypes])
    end
end


""" Import human or murine affinity data. """
@memoize function importKav(; murine = true, c1q = false, IgG2bFucose = false, retdf = false)
    if murine
        df = CSV.File(joinpath(dataDir, "murine-affinities.csv"), comment = "#") |> DataFrame!
    else
        df = CSV.File(joinpath(dataDir, "human-affinities.csv"), comment = "#") |> DataFrame!
    end

    IgGlist = copy(murine ? murineIgG : humanIgG)
    FcRecep = copy(murine ? murineFcgR : humanFcgR)
    if IgG2bFucose
        append!(IgGlist, [:IgG2bFucose])
    end
    if c1q
        append!(FcRecep, [:C1q])
    end
    df = stack(df; variable_name = :IgG, value_name = :Kav)
    df = unstack(df, :FcgR, :Kav)
    df[!, :IgG] = map(Symbol, df[!, :IgG])
    dropmissing!(df)
    df = df[in(IgGlist).(df.IgG), :]

    if retdf
        return df[!, [:IgG; FcRecep]]
    else
        return convert(Matrix{Float64}, df[!, FcRecep])
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

    df = CSV.File(joinpath(dataDir, filename), delim = ",", comment = "#") |> DataFrame!
    df[!, :Condition] .= Symbol.(df[!, :Condition])
    df[!, :Target] = 1.0 .- df[!, :Target] ./ 100.0
    if :Neutralization in propertynames(df)
        neut = -log.(df[!, :Neutralization] / 50.0)
        df[!, :Neutralization] .= replace!(neut, Inf => 0.0)
    end

    affinity = importKav(murine = true, c1q = c1q, IgG2bFucose = true, retdf = true)
    df = leftjoin(df, affinity, on = :Condition => :IgG)

    # The mG053 antibody doesn't bind to the virus
    if dataType == "HIV"
        df[df[:, :Label] .== "mG053", [:FcgRI, :FcgRIIB, :FcgRIII, :FcgRIV]] .= 0.0
    end

    df[df[:, :Background] .== "R1KO", :FcgRI] .= 0.0
    df[df[:, :Background] .== "R2KO", :FcgRIIB] .= 0.0
    df[df[:, :Background] .== "R3KO", :FcgRIII] .= 0.0
    df[df[:, :Background] .== "R1/3KO", [:FcgRI, :FcgRIII]] .= 0.0
    df[df[:, :Background] .== "R1/4KO", [:FcgRI, :FcgRIV]] .= 0.0
    df[df[:, :Background] .== "R4block", :FcgRIV] .= 0.0
    df[df[:, :Background] .== "gcKO", [:FcgRI, :FcgRIIB, :FcgRIII, :FcgRIV]] .= 0.0
    df[df[:, :Condition] .== :IgG1D265A, [:FcgRI, :FcgRIIB, :FcgRIII, :FcgRIV]] .= 0.0

    for pair in ["R" => "FcγR", "1" => "I", "2" => "II", "3" => "III", "4" => "IV", "gc" => "γc"]
        df[!, :Background] = map(x -> replace(x, pair), df.Background)
    end
    return df
end


""" Humanized mice data from Lux 2014, Schwab 2015 """
function importHumanized(dataType)
    if dataType in ["blood", "spleen", "bone"]
        df = CSV.File(joinpath(dataDir, "lux_humanized_CD19.csv"), delim = ",", comment = "#") |> DataFrame!
        df = dropmissing(df, Symbol(dataType), disallowmissing = true)
        df[!, :Target] = 1.0 .- df[!, Symbol(dataType)] ./ 100.0
        df[!, :Condition] .= :IgG1
        df = df[!, [:Genotype, :Concentration, :Condition, :Target]]
        affinity = importKav(murine = false, c1q = true, retdf = true)
    elseif dataType == "ITP"
        df = CSV.File(joinpath(dataDir, "schwab_ITP_humanized.csv"), delim = ",", comment = "#") |> DataFrame!
        df = stack(df, [:IgG1, :IgG2, :IgG3, :IgG4])
        df = disallowmissing!(df[completecases(df), :])
        rename!(df, [:variable => :Condition, :value => :Target])
        df[!, :Condition] .= Symbol.(df.Condition)

        df[!, :Target] .= 1.0 .- df.Target ./ 100.0
        df[!, :Donor] .= Symbol.(df.Donor)
        affinity = importKav(murine = false, c1q = false, retdf = true)
    else
        @error "Data type not found"
    end

    df = leftjoin(df, affinity, on = :Condition => :IgG)
    return df
end
