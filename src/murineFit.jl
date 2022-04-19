""" Import Apr 2022 murine in vitro data """
function importMurineInVitro(fn = "CHO-mFcgR-apr2022.csv")
    df = CSV.File(joinpath(dataDir, fn), comment = "#") |> DataFrame
    df[!, ["IgG1", "IgG2b", "IgG2c"]] .-= df[!, "TNP-BSA"]
    df = stack(df[!, Not("TNP-BSA")], ["IgG1", "IgG2b", "IgG2c"])
    rename!(df, "variable" => "Subclass")
    rename!(df, "value" => "Value")
    df = dropmissing(df[df."Receptor" .!= "CHO", :])
    return sort!(df, ["Receptor", "Subclass"])
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
    return polyfc(1e-9, KxStar, dfr."Valency", [recepExp[dfr."Receptor"] * dfr."Expression" / 100], 
        [1.0], reshape([dfr."Affinity"], 1, 1)).Lbound
end

function predictMurine(df::AbstractDataFrame; 
        Kav = importKav(; murine = true, retdf = true),
        conv = 1.62e4,
        kwargs...)
    # Add affinity information to df
    Kav[Kav."IgG" .== "IgG2a", "IgG"] .= "IgG2c"
    dft = innerjoin(df, Kav, on = "Subclass" => "IgG")
    dft = stack(dft, murineFcgR)
    rename!(dft, "value" => "Affinity")
    dft = dft[dft."Receptor" .== dft."variable", Not("variable")]
    sort!(dft, ["Receptor", "Subclass"])
    @assert df."Value" == dft."Value"   # must ensure the right order
    preds = Vector(undef, size(dft)[1])
    @Threads.threads for i = 1:size(dft)[1]
        preds[i] = predictMurine(dft[i, :]; kwargs...)
    end
    dft."Predict" = preds
    @assert all(isfinite(dft[!, "Predict"]))
    dft[!, "Predict"] ./= conv
    dft[dft."Predict" .< 1.0, "Predict"] .= 1.0
    return dft
end

@memoize function murineKavDist()
    Kav = importKav(; murine = true, retdf = true)
    function retDist(x)
        x = maximum([1e4, x])
        return inferLogNormal(x, x * 10)
    end
    Kav[!, Not("IgG")] = retDist.(Kav[!, Not("IgG")])
    return Kav
end

@model function murineFit(df, values)
    # Rtot sample
    Rtot = Vector(undef, length(murineFcgR))
    for ii in eachindex(Rtot)
        Rtot[ii] ~ truncated(LogNormal(18, 3), 100, 1E8)
        # Treat receptor amount as unknown with a wide prior
        # peak at ~10^5 with width ~10^4 to ~10^(6.5)
    end
    Rtotd = Dict([murineFcgR[ii] => Rtot[ii] for ii = 1:length(Rtot)])

    # Kav sample
    Kav_dist = Matrix(murineKavDist()[:, Not("IgG")])
    Kavd = importKav(; murine = true, retdf = true)
    Kav = Matrix(undef, size(Kav_dist)...)
    for ii in eachindex(Kav)
        Kav[ii] ~ truncated(Kav_dist[ii], 1e2, 1E10)
    end
    Kavd[!, Not("IgG")] = Kav

    KxStar ~ truncated(KxStarDist, 1E-18, 1E-9)
    # conversion factor: subtraction = 0.164
    conv ~ truncated(LogNormal(log(0.164), 2), 1e-6, 1e3)

    # fit predictions
    if all(Rtot .> 0.0) && all(Kav .> 0.0)
        df = predictMurine(deepcopy(df); recepExp = Rtotd, Kav = Kavd, KxStar = KxStar, conv = conv)
    else
        df = deepcopy(df)
        df."Predict" .= -1000.0
    end
    
    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end
