using DataFrames
using CSV
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
        df = CSV.read(joinpath(dataDir, "murine-FcgR-abundance.csv"))
    else
        df = CSV.read(joinpath(dataDir, "human-FcgR-abundance.csv"))
    end
    df = aggregate(df, [:Cells, :Receptor], geocmean)
    df = unstack(df, :Receptor, :Cells, :Count_geocmean)
    df[!, :Receptor] = map(Symbol, df[!, :Receptor])
    if murine
        df = df[in(murineFcgR).(df.Receptor), :]
    else
        df[df[:, :Receptor] .== :FcgRIIC, :Receptor] .= Symbol("FcgRIIC-13N")

        generic_type = Symbol.(["FcgRIIA", "FcgRIIB", "FcgRIIIA"])
        prefixes = ["FcgRIIA-131", "FcgRIIB-232", "FcgRIIIA-158"]
        options = [['H', 'R'], ['I', 'T'], ['V', 'F']]
        for i = 1:3
            rowidx = findfirst(df[:, :Receptor] .== generic_type[i])
            if genotype[i] == options[i][1]
                df[rowidx, :Receptor] = Symbol(prefixes[i] * options[i][1])
                insert!.(eachcol(df, false), rowidx + 1, [Symbol(prefixes[i] * options[i][2]); repeat([0.0], 7)])
            elseif genotype[i] == options[i][2]
                df[rowidx, :Receptor] = Symbol(prefixes[i] * options[i][2])
                insert!.(eachcol(df, false), rowidx, [Symbol(prefixes[i] * options[i][1]); repeat([0.0], 7)])
            else  # heterozygous
                insert!.(eachcol(df, false), rowidx, [Symbol(prefixes[i] * options[i][1]); Array(df[rowidx, 2:end]) ./ 2])
                insert!.(eachcol(df, false), rowidx + 1, [Symbol(prefixes[i] * options[i][2]); Array(df[rowidx, 2:end])])
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
        df = CSV.read(joinpath(dataDir, "murine-affinities.csv"), comment = "#")
    else
        df = CSV.read(joinpath(dataDir, "human-affinities.csv"), comment = "#")
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
    else
        @error "Data type not found"
    end

    df = CSV.read(joinpath(dataDir, filename), delim = ",", comment = "#")
    df[!, :Condition] .= Symbol.(df[!, :Condition])
    df[!, :Target] = 1.0 .- df[!, :Target] ./ 100.0
    if :Neutralization in names(df)
        neut = -log.(df[!, :Neutralization] / 50.0)
        df[!, :Neutralization] .= replace!(neut, Inf => 0.0)
    end

    affinity = importKav(murine = true, c1q = c1q, IgG2bFucose = true, retdf = true)
    df = join(df, affinity, on = :Condition => :IgG, kind = :left)

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
    return df
end


""" Humanized mice data from Lux 2014 """
function importHumanized(dataType)
    df = CSV.read(joinpath(dataDir, "lux_humanized_CD19.csv"), delim = ",", comment = "#")
    @assert dataType in ["blood", "spleen", "bone marrow"] "Data type not found"
    df = dropmissing(df, Symbol(dataType), disallowmissing = true)
    df[!, :Target] = 1.0 .- df[!, Symbol(dataType)] ./ 100.0
    df[!, :Condition] .= :IgG1
    df = df[!, [:Genotype, :Concentration, :Condition, :Target]]

    affinity = importKav(murine = false, c1q = true, retdf = true)
    df = join(df, affinity, on = :Condition => :IgG, kind = :left)
    return df
end



""" Import systems serology dataset. """
function importAlterMSG()
    dfF = CSV.read(joinpath(dataDir, "alter-MSB", "data-function.csv"))
    dfGP = CSV.read(joinpath(dataDir, "alter-MSB", "data-glycan-gp120.csv"))
    dfIGG = CSV.read(joinpath(dataDir, "alter-MSB", "data-luminex-igg.csv"))
    dfL = CSV.read(joinpath(dataDir, "alter-MSB", "data-luminex.csv"))
    dfMA = CSV.read(joinpath(dataDir, "alter-MSB", "meta-antigens.csv"))
    dfMD = CSV.read(joinpath(dataDir, "alter-MSB", "meta-detections.csv"))
    dfMG = CSV.read(joinpath(dataDir, "alter-MSB", "meta-glycans.csv"))
    dfMS = CSV.read(joinpath(dataDir, "alter-MSB", "meta-subjects.csv"))

    df = meltdf(dfL, view = true)
    newdfL = DataFrame(Rec = String[], Vir = String[], Sig = String[], Value = Float64[], Subject = Int64[])

    # Split column name into constituent parts
    for i = 1:size(df, 1)
        Ar = split(string(df.variable[i]), "."; limit = 3)
        if length(Ar) == 3
            push!(newdfL, [Ar[1], Ar[2], Ar[3], df.value[i], df.Column1[i]])
        else
            push!(newdfL, [Ar[1], Ar[2], "N/A", df.value[i], df.Column1[i]])
        end
    end

    return newdfL
end
