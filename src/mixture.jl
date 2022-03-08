using MultivariateStats
using StatsBase
import Statistics: cor

""" Load mixture in vitro binding data """
@memoize function loadMixData(fn = "lux_mixture_mar2021.csv")
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame

    df = stack(df, Not(["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])
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

lower(x) = quantile(x, 0.25)
upper(x) = quantile(x, 0.75)

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
    return combine(
        groupby(df, ["Valency", "Cell", "subclass_1", "subclass_2", "%_1", "%_2"]),
        valname => geocmean => valname,
        valname => geocstd => "std",
        valname => StatsBase.median => "Median",
        valname => lower => "xmin",
        valname => upper => "xmax",
    )
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
    return sort!(ndf, names(df)[in(["Valency", "Cell", "subclass_1", "subclass_2", "Experiment", "%_2"]).(names(df))])
end


""" 
General function to make subplots for every cell types and IgG pairs
splot() is a function that take dataframe with only a single cell type and IgG pair
    and output a plot
"""
function plotMixSubplots(splot::Function, df = loadMixData(); kwargs...)
    setGadflyTheme()

    cells = unique(df."Cell")
    pairs = unique(df[df."subclass_2" .!= "None", ["subclass_1", "subclass_2"]])
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

function R2(Actual, Predicted; logscale = true)
    if logscale
        return cor(log10.(Actual), log10.(Predicted))^2
    else
        return cor(Actual, Predicted)^2
    end
end


""" Four predictMix() below provide model predictions"""
function predictMix(
    cell,
    val::Real,
    IgGXname,
    IgGYname,
    IgGX,
    IgGY;
    recepExp = measuredRecepExp,
    KxStar = KxConst,
    Lbound = true,
    Kav::DataFrame = importKav(; murine = false, retdf = true),
    kwargs...,
)::Real
    IgGC = zeros(size(humanIgG))
    IgGC[IgGXname .== humanIgG] .= IgGX
    IgGC[IgGYname .== humanIgG] .= IgGY

    Kav = Matrix(Kav[!, [cell]])
    if IgGC' * Kav * [recepExp[cell]] <= 0.0
        return 0.0
    end
    res = try
        if Lbound
            polyfc(1e-9, KxStar, val, [recepExp[cell]], IgGC, Kav).Lbound
        else
            polyfc(1e-9, KxStar, val, [recepExp[cell]], IgGC, Kav).Rmulti
        end
    catch e
        println("Failed at predictMix():\n f = $val\n Rtot = $([recepExp[cell]])\n IgGC = $IgGC\n Kav = $Kav\n")
        rethrow(e)
    end
    return res
end

function predictMix(dfrow::DataFrameRow, IgGXname, IgGYname, IgGX, IgGY; kwargs...)
    val = "NewValency" in names(dfrow) ? dfrow."NewValency" : dfrow."Valency"
    return predictMix(dfrow."Cell", val, IgGXname, IgGYname, IgGX, IgGY; kwargs...)
end

predictMix(dfrow::DataFrameRow; kwargs...) = predictMix(dfrow, dfrow."subclass_1", dfrow."subclass_2", dfrow."%_1", dfrow."%_2"; kwargs...)

function predictMix(df::DataFrame; kwargs...)
    """ Will return another df object. """
    # Setup column
    df[!, "Predict"] .= predictMix(df[1, :]; kwargs...)
    for i = 2:size(df)[1]
        df[i, "Predict"] = predictMix(df[i, :]; kwargs...)
    end
    @assert all(isfinite(df[!, "Predict"]))
    df[df."Predict" .< 1.0, "Predict"] .= 1.0
    return df
end

""" PCA of isotype/combination x receptor matrix """
function mixtureDataPCA(; val = 0)
    df = averageMixData(loadMixData(); combSingle = true)
    if val > 0
        df = df[df."Valency" .== val, :]
    end
    id_cols = ["Valency", "subclass_1", "subclass_2", "%_1", "%_2"]
    wide = unstack(df, id_cols, "Cell", "Value")
    mat = Matrix(wide[!, Not(id_cols)])
    mat = coalesce.(mat, 0)
    M = MultivariateStats.fit(PCA, mat'; maxoutdim = 4)
    vars = principalvars(M)
    vars_expl = [sum(vars[1:i]) for i = 1:length(vars)] ./ tvar(M)

    score = MultivariateStats.transform(M, mat')'
    wide[!, "PC 1"] = score[:, 1]
    wide[!, "PC 2"] = score[:, 2]
    wide[!, "PC 3"] = score[:, 3]
    loading = projection(M)
    score_df = wide[!, vcat(id_cols, ["PC 1", "PC 2", "PC 3"])]
    loading_df = DataFrame("Cell" => unique(df."Cell"), "PC 1" => loading[:, 1], "PC 2" => loading[:, 2], "PC 3" => loading[:, 3])
    if "None" in df."subclass_2"
        score_df = combSing2pair(score_df)
    end
    score_df."Subclass Pair" = score_df."subclass_1" .* "-" .* score_df."subclass_2"
    return score_df, loading_df, vars_expl
end

#=
function mixtureANOVA()
    ## GLM and ANOVA are not included in the package. Results will be saved separately after run
    using GLM
    import ANOVA: anova
    df = loadMixData()
    df."Measurement" = string.(df."Valency") .* df."Cell" .* df."subclass_1" .* " " .* 
            string.(df."%_1") .* ", " .* df."subclass_2" .* " " .* string.(df."%_2")
    df."logValue" = log.(df."Value")

    model = fit(LinearModel,
            @formula(logValue ~ Measurement),
            df,
            contrasts = Dict(:Measurement => EffectsCoding()))
    return anova(model)
end
=#
