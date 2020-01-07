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
humanFcgR = [:FcgRI, :FcgRIIA, :FcgRIIB, :FcgRIIC, :FcgRIIIA, :FcgRIIIB]
murineActI = [1, -1, 1, 1]
humanActI = [1, 1, -1, 1, 1, 1]
dataDir = joinpath(dirname(pathof(FcgR)), "..", "data")

function humanFcgR_genotype(genotype)
    @assert length(genotype) == 3
    receps = Symbol.(["FcgRI", "FcgRIIC-13N", "FcgRIIIB"])
    append!(receps, Symbol.(["FcgRIIA-131" * genotype[1], "FcgRIIB-232" * genotype[2], "FcgRIIIA-158" * genotype[3]]))
    return sort(receps)
end

function importRtot(; murine = true)
    if murine
        df = CSV.read(joinpath(dataDir, "murine-FcgR-abundance.csv"))
    else
        df = CSV.read(joinpath(dataDir, "human-FcgR-abundance.csv"))
    end
    df = aggregate(df, [:Cells, :Receptor], geocmean)
    df = unstack(df, :Receptor, :Cells, :Count_geocmean)
    df[!, :Receptor] = map(Symbol, df[!, :Receptor])
    df = df[in(murine ? murineFcgR : humanFcgR).(df.Receptor), :]
    return convert(Matrix{Float64}, df[!, cellTypes])
end



""" Import human or murine affinity data. """
function importKav(; murine = true, c1q = false, retdf = false, IgG2bFucose = false, genotype = "RTF")
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
    df = df[!, murine ? murineFcgR : humanFcgR_genotype(genotype)]

    if retdf
        return df
    else
        return convert(Matrix{Float64}, df)
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
