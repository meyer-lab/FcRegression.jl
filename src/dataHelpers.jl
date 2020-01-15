using DataFrames
using CSV
import StatsBase.geomean

const KxConst = 6.31e-13 # 10^(-12.2)

function geocmean(x)
    x = convert(Vector, x)
    x[x .<= 1.0] .= 1.0
    return geomean(x)
end

cellTypes = [:ncMO, :cMO, :NKs, :Neu, :EO]
murineIgG = [:IgG1, :IgG2a, :IgG2b, :IgG3]
humanIgG = [:IgG1, :IgG2, :IgG3, :IgG4]
murineFcgR = [:FcgRI, :FcgRIIB, :FcgRIII, :FcgRIV]
humanFcgR =
    Symbol.(["FcgRI", "FcgRIIA-131H", "FcgRIIA-131R", "FcgRIIB-232I", "FcgRIIB-232T", "FcgRIIC-13N", "FcgRIIIA-158V", "FcgRIIIA-158F", "FcgRIIIB"])
murineActI = [1, -1, 1, 1]
humanActI = [1, 1, 1, -1, -1, 1, 1, 1, 1]
dataDir = joinpath(dirname(pathof(FcgR)), "..", "data")

function importRtot(; murine = true, genotype = "HIV", retdf = false)
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
function importKav(; murine = true, c1q = false, IgG2bFucose = false, retdf = false)
    if murine
        df = CSV.read(joinpath(dataDir, "murine-affinities.csv"), comment = "#")
    else
        df = CSV.read(joinpath(dataDir, "human-affinities.csv"), comment = "#")
    end

    if c1q == false
        df = filter(row -> row[:FcgR] != "C1q", df)
    end

    IgGlist = copy(murine ? murineIgG : humanIgG)
    if IgG2bFucose
        append!(IgGlist, [:IgG2bFucose])
    end
    df = stack(df; variable_name = :IgG, value_name = :Kav)
    df = unstack(df, :FcgR, :Kav)
    df = df[in(IgGlist).(df.IgG), :]

    if retdf
        return df
    else
        return convert(Matrix{Float64}, df[!, murine ? murineFcgR : humanFcgR])
    end
end


""" Import cell depletion data. """
function importDepletion(dataType; c1q = false)
    if dataType == "ITP"
        filename = "nimmerjahn-ITP.csv"
    elseif dataType == "blood"
        filename = "nimmerjahn-CD20-blood.csv"
    elseif dataType == "bone"
        filename = "nimmerjahn-CD20-bone.csv"
    elseif dataType == "melanoma"
        filename = "nimmerjahn-melanoma.csv"
    else
        @error "Data type not found"
    end

    df = CSV.read(joinpath(dataDir, filename), delim = ",", comment = "#")
    df[!, :Condition] = map(Symbol, df[!, :Condition])
    df[!, :Target] = 1.0 .- df[!, :Target] ./ 100.0

    affinityData = importKav(murine = true, c1q = c1q, retdf = true)
    df = join(df, affinityData, on = :Condition => :IgG, kind = :inner)

    df[df[:, :Background] .== "R1KO", :FcgRI] .= 0.0
    df[df[:, :Background] .== "R2KO", :FcgRIIB] .= 0.0
    df[df[:, :Background] .== "R3KO", :FcgRIII] .= 0.0
    df[df[:, :Background] .== "R1/3KO", [:FcgRI, :FcgRIII]] .= 0.0
    df[df[:, :Background] .== "R4block", :FcgRIV] .= 0.0
    df[df[:, :Background] .== "gcKO", [:FcgRI, :FcgRIIB, :FcgRIII, :FcgRIV]] .= 0.0
    return df
end
