using DataFrames
using Memoize
import CSV
using Distributions
using NLsolve
import Statistics: cor

const KxConst = 6.31e-13 # 10^(-12.2)

lower(x) = quantile(x, 0.25)
upper(x) = quantile(x, 0.75)

function geocmean(x)
    x = convert(Vector, x)
    x[x .<= 1.0] .= 1.0
    return geomean(x)
end


function R2(Actual, Predicted; logscale = true)
    if logscale
        return cor(log10.(Actual), log10.(Predicted))^2
    else
        return cor(Actual, Predicted)^2
    end
end

function bestfitline(xs, ys; logscale = false)
    if logscale
        xs = log10.(xs)
        ys = log10.(ys)
    end
    mx = mean(xs)
    my = mean(ys)
    dx = (xs .- mx)
    k = dx' * (ys .- my) / (dx' * dx)
    b = my - k * mx
    return k, b
end

const murineCellTypes = ["ncMO", "cMO", "NKs", "Neu", "EO", "Kupffer", "KupfferHi"]
const humanCellTypes = ["ncMO", "cMO", "NKs", "Neu", "EO"]
const murineIgG = ["IgG1", "IgG2a", "IgG2b", "IgG3"]
const murineIgGFucose = ["IgG1", "IgG2a", "IgG2b", "IgG3", "IgG2bFucose"]
const humanIgG = ["IgG1", "IgG2", "IgG3", "IgG4"]
const murineFcgR = ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"]
const humanFcgR = ["FcgRI", "FcgRIIA-131H", "FcgRIIA-131R", "FcgRIIB-232I", "FcgRIIIA-158F", "FcgRIIIA-158V", "FcgRIIIB"]
const humanFcgRiv = ["FcgRI", "FcgRIIA-131H", "FcgRIIA-131R", "FcgRIIB-232I", "FcgRIIIA-158F", "FcgRIIIA-158V"]

const murineActI = NamedArray([1.0, -1, 1, 1], ["FcgRI", "FcgRIIB", "FcgRIII", "FcgRIV"], "mFcgR")
const humanActI = NamedArray([1.0, 1, -1, 1, 1], ["FcgRI", "FcgRIIA", "FcgRIIB", "FcgRIIIA", "FcgRIIIB"], "hFcgR")


const colorReceptor = [
    colorant"hsl(195, 100%, 50%)",    # FcgRI
    colorant"hsl(56, 64%, 53%)",      # FcgRIIA-131H
    colorant"hsl(333, 100%, 71%)",    # FcgRIIA-131R
    colorant"hsl(166, 100%, 36%)",    # FcgRIIB-213I
    colorant"hsl(255, 100%, 83%)",    # FcgRIIIA-158F
    colorant"hsl(11, 100%, 58%)",     # FcgRIIIA-158V
    colorant"hsl(117, 90%, 68%)",     # FcgRIIIB
]
const colorSubclass = [
    colorant"hsl(115, 100%, 39%)",  # IgG1
    colorant"hsl(41, 100%, 39%)",   # IgG2
    colorant"hsl(352, 100%, 40%)",  # IgG3
    colorant"hsl(220, 100%, 45%)",  # IgG4
]
const colorValency = [
    colorant"#173f5f",    # 4
    colorant"#eb2200",    # 33
]
const colorAffinity = [
    colorant"navajowhite2",   # documented
    colorant"firebrick4",     # updated
]

const dataDir = joinpath(dirname(pathof(FcRegression)), "..", "data")

@memoize function importRtot_readcsv(; murine::Bool, genotype = "HIV", retdf = true, cellTypes::Union{Nothing, Vector, NamedVector} = nothing)
    if murine
        df = CSV.File(joinpath(dataDir, "murine-leukocyte-FcgR-abundance.csv"), comment = "#") |> DataFrame
    else
        df = CSV.File(joinpath(dataDir, "human-leukocyte-FcgR-abundance.csv"), comment = "#") |> DataFrame
    end
    if cellTypes === nothing
        cellTypes = murine ? murineCellTypes : humanCellTypes
    elseif cellTypes isa NamedArray
        cellTypes = names(cellTypes)[1]
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
    df = df[in(murine ? murineFcgR : humanFcgR).(df."Receptor"), :]
    @assert df.Receptor == (murine ? murineFcgR : humanFcgR)
    if retdf
        return df[!, ["Receptor"; names(df)[in(cellTypes).(names(df))]]]
    else
        return Matrix{Float64}(df[!, cellTypes])
    end
end

importRtot(; kwargs...) = deepcopy(importRtot_readcsv(; kwargs...))

""" Import human or murine affinity data. """
@memoize function importKav_readcsv(; murine::Bool, c1q = false, IgG2bFucose = false, retdf = true)
    if murine
        df = CSV.File(joinpath(dataDir, "murine-affinities.csv"), comment = "#") |> DataFrame
    else
        df = CSV.File(joinpath(dataDir, "human-affinities.csv"), comment = "#") |> DataFrame
    end

    IgGlist = copy(murine ? murineIgG : humanIgG)
    FcRecep = copy(murine ? murineFcgR : humanFcgR)
    if IgG2bFucose
        append!(IgGlist, ["IgG2bFucose", "IgG1SA", "IgG2bSA"])
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

importKav(; kwargs...) = deepcopy(importKav_readcsv(; kwargs...))

""" A more accurate way to infer logNormal distribution with exact mode and IQR """
@memoize function inferLogNormal(mode, iqr)
    function logNormalParams!(f, v)
        f[1] = exp(v[1] - v[2]^2) - mode
        f[2] = 2 * exp(v[1]) * sinh(0.6745 * v[2]) - iqr
        # from Wikipedia and Dewey Lonzo Whaley (ETSU)'s thesis, Eq. 36
    end
    xs = nlsolve(logNormalParams!, [log(mode), 1.0]).zero
    return LogNormal(xs[1], xs[2])
end

""" Import measurements of receptor amounts. """
@memoize function importRtotDist_readcsv(dat::Symbol; regular = false, retdf = true)
    @assert dat in [:hCHO, :hRob, :mCHO, :mLeuk]
    if dat in [:hCHO, :hRob]
        df = if dat == :hRob
            CSV.File(joinpath(dataDir, "robinett_FcgR_quant.csv"), delim = ",", comment = "#") |> DataFrame
        else
            CSV.File(joinpath(dataDir, "CHO_FcgR_quant.csv"), delim = ",", comment = "#") |> DataFrame
        end
        res = if regular
            [geomean(df[df."Receptor" .== rcp, "Measurements"]) for rcp in unique(df."Receptor")]
        else
            [fit_mle(LogNormal, df[df."Receptor" .== rcp, "Measurements"]) for rcp in unique(df."Receptor")]
        end
        if retdf
            return Dict([humanFcgRiv[i] => res[i] for i = 1:length(res)])
        else
            return res
        end
    end
    if dat == :mCHO
        ref = 1e6
        res = [regular ? ref : inferLogNormal(ref, ref * 10) for ii = 1:length(murineFcgR)]
        if retdf
            return Dict([murineFcgR[i] => res[i] for i = 1:length(res)])
        else
            return res
        end
    elseif dat == :mLeuk
        cells = unique(importMurineLeukocyte(; average = false)."ImCell")
        df = CSV.File(joinpath(dataDir, "murine-FcgR-abundance.csv"), comment = "#") |> DataFrame
        df[df."Count" .< 1.0, "Count"] .= 1.0

        function findDist(x)
            if regular
                return geomean(x)
            end
            if length(x) <= 1
                return inferLogNormal(x[1], x[1])
            end
            return inferLogNormal(geomean(x), maximum([quantile(x, 0.7) - quantile(x, 0.3), 5.0]))
        end

        ndf = combine(groupby(df, ["Cells", "Receptor"]), "Count" => findDist => "Distribution")
        ndf = dropmissing(unstack(ndf, "Receptor", "Cells", "Distribution"))
        @assert ndf."Receptor" == murineFcgR
        ndf = ndf[!, ["Receptor"; names(ndf)[in(cells).(names(ndf))]]]
        if retdf
            return ndf
        else
            return Matrix(ndf[!, Not("Receptor")])
        end
    end
end

importRtotDist(dat; kwargs...) = deepcopy(importRtotDist_readcsv(dat; kwargs...))

@memoize function importKavDist_readcsv(; murine::Bool, regularKav = false, retdf = true, CD16b = false)
    local Kav
    if murine
        Kav = importKav(; murine = true, retdf = true)
        Kav = Kav[Kav."IgG" .!= "IgG3", :]
        function retDist(x; regularKav = regularKav)
            x = maximum([1e4, x])
            if regularKav
                return x
            end
            return inferLogNormal(x, x)
        end
        Kav[Kav."IgG" .== "IgG2a", "IgG"] .= "IgG2c"
        Kav[!, Not("IgG")] = retDist.(Kav[!, Not("IgG")], regularKav = regularKav)
        sort!(Kav, "IgG")
    else # human
        df = CSV.File(joinpath(dataDir, "human_affinities_variance.csv"), delim = ",", comment = "#") |> DataFrame
        function parstr(x, regularKav = false)
            params = parse.(Float64, split(x, "|"))
            params .*= 1e5      # Bruhns data is written in 1e5 units
            if regularKav
                return params[1]
            end
            params[1] = maximum([params[1], 1e4])   # minimum affinity as 1e4 M-1
            params[2] = maximum([params[2], 1e5])   # minimum IQR as 1e5 M-1
            return inferLogNormal(params[1], params[2])
        end
        Kav = parstr.(df[:, Not("IgG")], regularKav)
        insertcols!(Kav, 1, "IgG" => df[:, "IgG"])
        if any(startswith.(names(Kav), "FcgRIIIB")) && !CD16b
            Kav = Kav[!, Not(startswith.(names(Kav), "FcgRIIIB"))]
        end
    end
    if retdf
        return Kav
    else
        return Matrix(Kav[:, Not("IgG")])
    end
end

importKavDist(; kwargs...) = deepcopy(importKavDist_readcsv(; kwargs...))

""" Humanized mice data from Lux 2014, Schwab 2015 """
function importHumanized()
    df = CSV.File(joinpath(dataDir, "humanized_mice_ITP.csv"), delim = ",", comment = "#") |> DataFrame
    df = stack(df, ["IgG1", "IgG2", "IgG3", "IgG4"])
    df = disallowmissing!(df[completecases(df), :])
    rename!(df, ["variable" => "Condition", "value" => "Target"])
    df[!, "Target"] .= 1.0 .- df.Target ./ 100.0
    df."Genotype" = [g[1] * "I" * g[3] for g in df."Genotype"]   # FcgRIIB default as 232I
    @assert maximum(df."Target") <= 1.0
    return df
end


function printDocumentedKavSupTable()
    df = CSV.File(joinpath(dataDir, "FcgR-Ka-Bruhns_with_variance.csv"), delim = ",", comment = "#") |> DataFrame
    function parstr(x)
        params = parse.(Float64, split(x, "|"))
        params .*= 1e5      # Bruhns data is written in 1e5 units
        params[1] = maximum([params[1], 1e4])   # minimum affinity as 1e4 M-1
        params[2] = maximum([params[2], 1e5])   # minimum variance as 1e5 M-1
        dist = inferLogNormal(params[1], params[2])
        return @sprintf("%.2e ± %.2e\n~lnN(μ=%.2f, σ=%.2f)", params[1], params[2], dist.μ, dist.σ)
    end
    CSV.write("documented_Kav.csv", parstr.(df[:, Not("IgG")]))
end

function printFittedKavSupTable()
    c = rungMCMC("humanKavfit_0701.dat"; dat = :hCHO, mcmc_iter = 1_000)
    pnames = [String(s) for s in c.name_map[1]]

    function parchain(s)
        return @sprintf("%.3e\n%.3e~%.3e", median(c[s].data), quantile(c[s].data[:], 0.25), quantile(c[s].data[:], 0.75))
    end

    Kavd = importKavDist(; murine = false, regularKav = true, retdf = true)
    Kav = [parchain("Kav[$i]") for i = 1:sum(startswith.(pnames, "Kav"))]
    Kavd[!, Not("IgG")] = typeof(Kav[1, 1]).(reshape(Kav, size(Kavd)[1], :))
    CSV.write("updated_Kav.csv", Kavd)
end
