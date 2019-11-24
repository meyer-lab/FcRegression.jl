"""
Guidance on DataFrame handling:
1. First thing first: DataFrame is different from Matrix!
2. All saved .csv files must be a valid dataframe (explained below) when
    imported by a one-line read command (delim = ',', comment = "#", etc.)
3. No two columns shall bear the same name (especially in csv files)
4. Keep the row indices (1,2,3,...) auto-generated by the import function
5. Use long format to store those with replicates
6. Use consistent name in all .csv files
    (e.g. never have both "FcgRI" and "FcgR1", even in different files)
7. Check/sort/reorder the order of rows and cols whenever that matters
"""

using DataFrames
using CSV

function geocmean(x)
    x = convert(Vector, x)
    @. x[x <= 1.0] = 1.0
    return exp( sum(log.(x))/length(x) )
end

cellTypes = [:ncMO, :cMO, :NKs, :Neu, :EO]

function Rpho_mouse()
    df = CSV.read("../data/murine-FcgR-abundance.csv")
    df = aggregate(df, [:Cells, :Receptor], geocmean)
    df = unstack(df, :Receptor, :Cells, :Count_geocmean)
    return convert(Matrix{Float64}, df[!, cellTypes])
end


""" Import human or murine affinity data. """
function importKav(; murine=true, c1q=false)
    if murine
        df = CSV.read("../data/murine-affinities.csv", comment="#")
    else
        df = CSV.read("../data/human-affinities.csv", comment="#")
    end

    if c1q == false
        df = filter(row -> row[:FcgR] != "C1q", df)
    end

    df = melt(df; variable_name=:IgG, value_name=:Kav)
    df = unstack(df, :FcgR, :Kav)

    return df
end


""" Import cell depletion data. """
function import_depletion(dataType; c1q=false)
    if dataType == "melanoma"
        df = CSV.read("../data/nimmerjahn-melanoma.csv", comment="#")
    end

    df[!, :Condition] = map(Symbol, df[!, :Condition])

    affinityData = importKav(murine=true, c1q=c1q)

    df = join(df, affinityData, on = :Condition => :IgG, kind = :outer)

    df[df[:, :Background] .== "R1KO", :FcgRI] .= 0.0
    df[df[:, :Background] .== "R2KO", :FcgRIIB] .= 0.0
    df[df[:, :Background] .== "R3KO", :FcgRIII] .= 0.0
    df[df[:, :Background] .== "R1/3KO", [:FcgRI, :FcgRIII]] .= 0.0
    df[df[:, :Background] .== "R4block", :FcgRIV] .= 0.0
    df[df[:, :Background] .== "gcKO", [:FcgRI, :FcgRIIB, :FcgRIII, :FcgRIV]] .= 0.0

    return df
end


const KxConst = 6.31e-13 # 10^(-12.2)
