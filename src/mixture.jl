using MultivariateStats
using Impute
using StatsBase
import Statistics: cor

""" Load mixture in vitro binding data """
function loadMixData(fn = "lux_mixture_mar2021.csv"; discard_small = false)
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame

    df = stack(df, Not(["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])
    if discard_small
        df = df[df[!, "Value"] .> 100, :]   # discard small measurements
    end
    df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0

    df[!, "%_1"] ./= 100.0
    df[!, "%_2"] ./= 100.0

    replace!(df."Cell", "CHO-hFcgRIIA-131His" => "FcgRIIA-131H")
    replace!(df."Cell", "CHO-hFcgRIIB" => "FcgRIIB-232I")
    replace!(df."Cell", "CHO-hFcgRIIIA-131Val" => "FcgRIIIA-158V")
    replace!(df."Cell", "CHO-FcgRIA" => "FcgRI")
    replace!(df."Cell", "CHO-hFcgRIIA-131Arg" => "FcgRIIA-131R")
    replace!(df."Cell", "CHO-hFcgRIIIA-158Phe" => "FcgRIIIA-158F")

    return sort!(df, ["Valency", "Cell", "subclass_1", "subclass_2", "Experiment", "%_2"])
end

""" Make statistics of individual cell types and subclass types """
function averageMixData(df = loadMixData())
    lower(x) = quantile(x, 0.25)
    upper(x) = quantile(x, 0.75)
    valname = "Adjusted" in names(df) ? "Adjusted" : "Value"
    return combine(
        groupby(df, ["Valency", "Cell", "subclass_1", "subclass_2", "%_1", "%_2"]),
        valname => geocmean => valname,
        valname => geocstd => "std",
        valname => StatsBase.median => "Median",
        valname => lower => "xmin",
        valname => upper => "xmax",
    )
end

""" 
General function to make subplots for every cell types and IgG pairs
splot() is a function that take dataframe with only a single cell type and IgG pair
    and output a plot
"""
function plotMixSubplots(splot::Function, df = loadMixData(); avg = false, kwargs...)
    setGadflyTheme()
    if avg
        df = averageMixData(df)
    end

    cells = unique(df."Cell")
    pairs = unique(df[!, ["subclass_1", "subclass_2"]])
    lcells = length(cells)
    lpairs = size(pairs, 1)
    pls = Vector(undef, lcells * lpairs)
    palette = [Scale.color_discrete().f(3)[1], Scale.color_discrete().f(3)[3], Scale.color_discrete().f(3)[2], Scale.color_discrete().f(3)[4:end]...]

    for (i, pairrow) in enumerate(eachrow(pairs))
        for (j, cell) in enumerate(cells)
            IgGXname, IgGYname = pairrow."subclass_1", pairrow."subclass_2"
            ndf = df[(df."Cell" .== cell) .& (df."subclass_1" .== IgGXname) .& (df."subclass_2" .== IgGYname), :]
            pls[(j - 1) * lpairs + (i - 1) + 1] = splot(ndf; kwargs...)
        end
    end
    return plotGrid((lcells, lpairs), pls; sublabels = false)

end


const measuredRecepExp = Dict(
    "FcgRI" => 101493.689,
    "FcgRIIA-131H" => 1006302.484,
    "FcgRIIA-131R" => 190432.6753,
    "FcgRIIB-232I" => 75085.07599,
    "FcgRIIIA-158F" => 634324.0675,
    "FcgRIIIA-158V" => 979451.9884,
)  # geometric mean precalculated


function ols(Actual, Predicted; logscale = true)
    if logscale
        return log.(Actual) \ log.(Predicted)
    end
    return Actual \ Predicted
end

function R2(Actual, Predicted)
    return cor(log10.(Actual), log10.(Predicted))^2.0
end

""" Three predictMix() below provide model predictions"""
function predictMix(dfrow::DataFrameRow, IgGXname, IgGYname, IgGX, IgGY; recepExp = measuredRecepExp, KxStar = KxConst)
    IgGC = zeros(size(humanIgG))
    IgGC[IgGXname .== humanIgG] .= IgGX
    IgGC[IgGYname .== humanIgG] .= IgGY

    Kav = importKav(; murine = false, retdf = true)
    Kav = Matrix(Kav[!, [dfrow."Cell"]])
    val = "NewValency" in names(dfrow) ? dfrow."NewValency" : dfrow."Valency"
    res = try
        polyfc(1e-9, KxStar, val, [recepExp[dfrow."Cell"]], IgGC, Kav).Lbound
    catch e
        println(val, [recepExp[dfrow."Cell"]], IgGC, Kav)
        rethrow(e)
    end
    return res
end

predictMix(dfrow::DataFrameRow; recepExp = measuredRecepExp, KxStar = KxConst) =
    predictMix(dfrow, dfrow."subclass_1", dfrow."subclass_2", dfrow."%_1", dfrow."%_2"; recepExp = recepExp, KxStar = KxStar)

function predictMix(df::DataFrame; recepExp = measuredRecepExp, KxStar = KxConst)
    """ will return another df object """
    df = copy(df)
    df[!, "Predict"] .= convert(typeof(KxStar), 1.0)
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictMix(df[i, :]; recepExp = recepExp, KxStar = KxStar)
    end
    df[df."Predict" .< 1.0, "Predict"] .= 1.0
    return df
end

""" PCA of isotype/combination x receptor matrix """
function mixtureDataPCA(; val = 0)
    df = averageMixData(loadMixData(; discard_small = false))
    if df > 0
        df = df[df."Valency" .== val, :]
    end
    id_cols = ["Valency", "subclass_1", "subclass_2", "%_1", "%_2"]
    wide = unstack(df, id_cols, "Cell", "Value")
    mat = Matrix(wide[!, Not(id_cols)])
    mat = coalesce.(mat, 0)
    M = fit(PCA, mat'; maxoutdim = 4)
    vars = principalvars(M)
    vars_expl = [sum(vars[1:i]) for i = 1:length(vars)] ./ tvar(M)

    score = MultivariateStats.transform(M, mat')'
    wide[!, "PC 1"] = score[:, 1]
    wide[!, "PC 2"] = score[:, 2]
    loading = projection(M)
    score_df = wide[!, vcat(id_cols, ["PC 1", "PC 2"])]
    loading_df = DataFrame("Cell" => unique(df."Cell"), "PC 1" => loading[:, 1], "PC 2" => loading[:, 2])
    return score_df, loading_df, vars_expl
end
