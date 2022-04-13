function importMurineInVitro(fn = "CHO-mFcgR-apr2022.csv"; subTdivF = true)
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame
    df = stack(df, ["TNP-BSA", "IgG1", "IgG2b", "IgG2c"])
    rename!(df, "variable" => "Subclass")
    rename!(df, "value" => "Value")
    dfCHO = dropmissing(df[df."Receptor" .== "CHO", Not("Expression")])
    dfCHO = combine(groupby(dfCHO, "Subclass"), "Value" => geocmean => "BaseValue")
    dfn = dropmissing(df[df."Receptor" .!= "CHO", :])
    dfn = innerjoin(dfCHO, dfn, on = "Subclass")
    if subTdivF
        dfn."Value" = dfn."Value" .- dfn."BaseValue"
    else
        dfn."Value" = dfn."Value" ./ dfn."BaseValue"
    end
    return dfn[!, Not("BaseValue")]
end

function predictMurine(dfr::DataFrameRow; 
        KxStar = KxConst, 
        recepExp = Dict(
            "FcgRI" => 1e5,
            "FcgRIIB" => 1e5,
            "FcgRIII" => 1e5,
            "FcgRIV" => 1e5,
        ))
    if dfr."Subclass" == "TNP-BSA" || dfr."Affinity" <= 0.0
        return 0.0
    end
    return polyfc(1e-9, KxStar, dfr."Valency", [recepExp[dfr."Receptor"] * dfr."Expression"], 
        [1.0], reshape([dfr."Affinity"], 1, 1)).Lbound
end

function predictMurine(df::AbstractDataFrame; kwargs...)
    # Add affinity information to df
    Kav = importKav(;murine = true, retdf = true)
    Kav[Kav."IgG" .== "IgG2a", "IgG"] .= "IgG2c"
    push!(Kav, ["TNP-BSA", 0.0, 0.0, 0.0, 0.0])
    dft = innerjoin(df, Kav, on = "Subclass" => "IgG")
    dft = stack(dft, murineFcgR)
    rename!(dft, "value" => "Affinity")
    df = dft[dft."Receptor" .== dft."variable", Not("variable")]
    df."Predict" .= 0.0
    for i = 1:size(df)[1]
        df[i, "Predict"] = predictMurine(df[i, :]; kwargs...)
    end
    @assert all(isfinite(df[!, "Predict"]))
    df[df."Predict" .< 1.0, "Predict"] .= 1.0
    return df
end