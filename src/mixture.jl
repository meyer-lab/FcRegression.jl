using MultivariateStats
using Impute
using StatsBase
import Statistics: cor
using GLM

""" Load mixture in vitro binding data """
function loadMixData(fn = "lux_mixture_mar2021.csv";)
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame

    df = stack(df, Not(["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])
    df = df[df[!, "Value"] .> 100, :]   # discard small measurements
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
        valname => lower => "ymin",
        valname => upper => "ymax",
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
    return plotGrid((lcells, lpairs), pls; sublabel = false)

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
        df = DataFrame(A = log.(Actual), B = log.(Predicted))
    else
        df = DataFrame(A = Actual, B = Predicted)
    end
    return coef(lm(@formula(B ~ A + 0), df))[1]
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


function MixtureFitLoss(df, ValConv::Vector, ExpConv::Vector; logscale = false)
    """ Loss function to be optimized given the conversion factors, called only by MixtureFit() """
    df = copy(df)
    dtype = promote_type(eltype(df."Value"), eltype(ValConv), eltype(ExpConv))
    df[!, "Adjusted"] .= convert.(dtype, df[!, "Value"])
    df[df[!, "Adjusted"] .<= 1.0, "Adjusted"] .= 1.0

    @assert all(df[!, "Adjusted"] .>= 1.0)
    @assert all(df[!, "Predict"] .>= 1.0)
    if any(ValConv .<= 0.0) || any(ExpConv .<= 0.0)
        return Inf, df
    end

    loss = 0.0
    for (iv, val) in enumerate(unique(df."Valency"))
        for (ie, exp) in enumerate(unique(df."Experiment"))
            sdf = @view df[(df."Valency" .== val) .& (df."Experiment" .== exp), :]
            sdf."Adjusted" .*= ValConv[iv] * ExpConv[ie]
            if logscale
                loss += norm(log.(sdf."Adjusted") .- log.(sdf."Predict"), 2) / nrow(sdf)
            else
                loss += norm(sdf."Adjusted" .- sdf."Predict", 2) / nrow(sdf)
            end
        end
    end
    return loss, df
end


function MixtureFit(df; logscale = false)
    """ Two-way fitting for valency and experiment (day) """
    if !("Predict" in names(df))
        df = predictMix(df)
    end
    nv, ne = length(unique(df."Valency")), length(unique(df."Experiment"))
    f(p::Vector, q::Vector) = MixtureFitLoss(df, [1.0; p], q; logscale = logscale)[1]
    f(v::Vector) = f(v[1:(nv - 1)], v[nv:end])
    init_v = ones(nv + ne - 1)
    before_loss = MixtureFitLoss(df, [1.0; init_v[1:(nv - 1)]], init_v[nv:end]; logscale = logscale)[1]
    od = OnceDifferentiable(f, init_v; autodiff = :forward)
    res = optimize(od, init_v, BFGS()).minimizer
    p, q = [1.0; res[1:(nv - 1)]], res[nv:end]
    res = MixtureFitLoss(df, p, q; logscale = logscale)
    return Dict(
        "loss" => res[1],
        "df" => res[2],
        "ValConv" => Dict([(name, p[i]) for (i, name) in enumerate(unique(df."Valency"))]),
        "ExpConv" => Dict([(name, q[i]) for (i, name) in enumerate(unique(df."Experiment"))]),
    )
end
