using MultivariateStats
using StatsBase
import Statistics: cor

""" Load mixture in vitro binding data """
function loadMixData(fn = "lux_mixture_mar2021.csv")
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame

    df = stack(df, Not(["Valency", "Cell", "subclass_1", "%_1", "subclass_2", "%_2"]), variable_name = "Experiment", value_name = "Value")
    df = dropmissing(df)
    df[!, "Value"] = convert.(Float64, df[!, "Value"])
    df[(df[!, "Value"]) .< 1.0, "Value"] .= 1.0

    baseline = combine(groupby(df, "Experiment"), "Value" => geomean => "Baseline")
    df = innerjoin(df, baseline, on = "Experiment")
    df[!, "Value"] ./= df[!, "Baseline"]    # normalize fluorescence by daily geomean
    df = df[!, Not(["Experiment", "Baseline"])]
    rename!(df, "Cell" => "Receptor")

    df[!, "%_1"] ./= 100.0
    df[!, "%_2"] ./= 100.0

    replace!(df."Receptor", "CHO-hFcgRIIA-131His" => "FcgRIIA-131H")
    replace!(df."Receptor", "CHO-hFcgRIIB" => "FcgRIIB-232I")
    replace!(df."Receptor", "CHO-hFcgRIIIA-131Val" => "FcgRIIIA-158V")
    replace!(df."Receptor", "CHO-FcgRIA" => "FcgRI")
    replace!(df."Receptor", "CHO-hFcgRIIA-131Arg" => "FcgRIIA-131R")
    replace!(df."Receptor", "CHO-hFcgRIIIA-158Phe" => "FcgRIIIA-158F")

    return sort!(df, ["Valency", "Receptor", "subclass_1", "subclass_2", "%_2"])
end

function importRobinett()
    df = CSV.File(joinpath(dataDir, "robinett/Luxetal2013-Fig2Bmod.csv"), delim = ",", comment = "#") |> DataFrame
    for i = 1:4
        cn = "Replicate $i"
        df[!, cn] ./= geomean(df[Not(ismissing.(df[!, cn])), cn])
    end
    df = dropmissing(stack(df, Not(["Receptor", "Antibody", "Valency"])))
    rename!(df, ["variable" => "Experiment", "value" => "Value"])
    rename!(df, ["Antibody" => "Subclass"])

    return sort!(df, ["Valency", "Receptor", "Subclass", "Experiment"])
end

function importMurineLeukocyte(fn = "leukocyte-apr2022.csv"; average = true)
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame
    #df = df[df."ImCell" .!= "Tcell", :]   # throw away T cells
    df = df[df."Experiment" .!= 7, :]   # throw away Exp. 7
    df[!, ["IgG1", "IgG2b", "IgG2c"]] .-= df[!, "TNP-BSA"]  # subtract TNP-BSA control from meas
    df = df[!, Not("TNP-BSA")]
    df = dropmissing(stack(df, ["IgG1", "IgG2b", "IgG2c"], variable_name = "Subclass", value_name = "Value"))
    df[df."Value" .< 1.0, "Value"] .= 1.0   # clip values to 1.0
    baseline = combine(groupby(df, "Experiment"), "Value" => geomean => "Baseline")
    df = innerjoin(df, baseline, on = "Experiment")
    df[!, "Value"] ./= df[!, "Baseline"]    # normalize fluorescence by daily geomean
    df = df[!, Not(["Experiment", "Baseline"])]
    if average
        df = combine(
            groupby(df, Not("Value")),
            "Value" => geomean => "Value",
            "Value" => (xs -> quantile(xs, 0.25)) => "xmin",
            "Value" => (xs -> quantile(xs, 0.75)) => "xmax",
        )
    end
    return sort!(df, ["ImCell", "Subclass", "Valency"])
end

""" Import Apr 2022 murine in vitro data """
function importMurineInVitro(fn = "CHO-mFcgR-apr2022.csv")
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame
    df = df[df."Experiment" .!= 5, :]   # throw away Exp. 5
    df[!, ["IgG1", "IgG2b", "IgG2c"]] .-= df[!, "TNP-BSA"]  # subtract TNP-BSA control from meas
    df = df[df."Receptor" .!= "CHO", Not("TNP-BSA")]
    df = dropmissing(stack(df, ["IgG1", "IgG2b", "IgG2c"], variable_name = "Subclass", value_name = "Value"))
    baseline = combine(groupby(df, "Experiment"), "Value" => geomean => "Baseline")
    df = innerjoin(df, baseline, on = "Experiment")
    df[!, "Value"] ./= df[!, "Baseline"]    # normalize fluorescence by daily geomean
    df = df[!, Not(["Experiment", "Baseline"])]
    return sort!(df, ["Receptor", "Subclass"])
end


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
    combining = ["Value"]
    if "Experiment" in names(df)
        append!(combining, ["Experiment"])
    end
    return combine(
        groupby(df, Not(combining)),
        valname => geomean => valname,
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
    return sort!(ndf, names(df)[in(["Valency", "Receptor", "subclass_1", "subclass_2", "Experiment", "%_2"]).(names(df))])
end


""" 
General function to make subplots for every cell types and IgG pairs
splot() is a function that take dataframe with only a single cell type and IgG pair
    and output a plot
"""
function plotMixSubplots(splot::Function, df = loadMixData(); kwargs...)
    setGadflyTheme()

    cells = unique(df."Receptor")
    pairs = unique(df[df."subclass_2" .!= "None", ["subclass_1", "subclass_2"]])
    lcells = length(cells)
    lpairs = size(pairs, 1)
    pls = Vector(undef, lcells * lpairs)

    for (i, pairrow) in enumerate(eachrow(pairs))
        for (j, cell) in enumerate(cells)
            IgGXname, IgGYname = pairrow."subclass_1", pairrow."subclass_2"
            ndf = df[(df."Receptor" .== cell) .& (df."subclass_1" .== IgGXname) .& (df."subclass_2" .== IgGYname), :]
            pls[(j - 1) * lpairs + (i - 1) + 1] = splot(ndf; kwargs...)
        end
    end
    return plotGrid((lcells, lpairs), pls; sublabels = false)

end


function R2(Actual, Predicted; logscale = true)
    if logscale
        return cor(log10.(Actual), log10.(Predicted))^2
    else
        return cor(Actual, Predicted)^2
    end
end

""" PCA of isotype/combination x receptor matrix """
function mixtureDataPCA(; val = 0)
    df = averageMixData(loadMixData(); combSingle = true)
    if val > 0
        df = df[df."Valency" .== val, :]
    end
    id_cols = ["Valency", "subclass_1", "subclass_2", "%_1", "%_2"]
    wide = unstack(df, id_cols, "Receptor", "Value")
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
    loading_df = DataFrame("Receptor" => unique(df."Receptor"), "PC 1" => loading[:, 1], "PC 2" => loading[:, 2], "PC 3" => loading[:, 3])
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
    df."Measurement" = string.(df."Valency") .* df."Receptor" .* df."subclass_1" .* " " .* 
            string.(df."%_1") .* ", " .* df."subclass_2" .* " " .* string.(df."%_2")
    df."logValue" = log.(df."Value")

    model = fit(LinearModel,
            @formula(logValue ~ Measurement),
            df,
            contrasts = Dict(:Measurement => EffectsCoding()))
    return anova(model)
end
=#
