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

""" Priors for murine affinities. Not including IgG3"""
@memoize function murineKavDist_nested(; regularKav = false, retdf = true)
    Kav = importKav(; murine = true, retdf = true)
    Kav = Kav[Kav."IgG" .!= "IgG3", :]
    function retDist(x; regularKav = regularKav)
        x = maximum([1e4, x])
        if regularKav
            return x
        end
        return inferLogNormal(x, x * 10)
    end
    Kav[Kav."IgG" .== "IgG2a", "IgG"] .= "IgG2c"
    Kav[!, Not("IgG")] = retDist.(Kav[!, Not("IgG")], regularKav = regularKav)
    sort!(Kav, "IgG")
    if !retdf
        return Matrix(Kav[!, Not("IgG")])
    end
    return Kav
end

murineKavDist(; kwargs...) = deepcopy(murineKavDist_nested(; kwargs...))

const InVitroMurineRcpExp = Dict("FcgRI" => 1e6, "FcgRIIB" => 1e6, "FcgRIII" => 1e6, "FcgRIV" => 1e6)

function predictMurine(dfr::DataFrameRow; KxStar = KxConst, recepExp = InVitroMurineRcpExp, f33 = 33)
    if (dfr."Subclass" == "TNP-BSA") || (dfr."Affinity" <= 0.0)
        return 0.0
    end
    if dfr."Affinity" >= Inf
        return Inf
    end
    @assert dfr."Valency" == 33     # use valency fitted from human; only deal with f = 33 here
    return polyfc(1e-9, KxStar, f33, [recepExp[dfr."Receptor"] * dfr."Expression" / 100], [1.0], reshape([dfr."Affinity"], 1, 1)).Lbound
end

function predictMurine(df::AbstractDataFrame; Kav = murineKavDist(; regularKav = true), kwargs...)
    # Add affinity information to df
    Kav = deepcopy(Kav)
    Kav[Kav."IgG" .== "IgG2a", "IgG"] .= "IgG2c"
    dft = innerjoin(df, Kav, on = "Subclass" => "IgG")
    dft = stack(dft, murineFcgR)
    rename!(dft, "value" => "Affinity")
    dft = dft[dft."Receptor" .== dft."variable", Not("variable")]
    sort!(dft, ["Receptor", "Subclass"])
    @assert df."Value" == dft."Value"   # must ensure the right order
    preds = Vector(undef, size(dft)[1])
    Threads.@threads for i = 1:size(dft)[1]
        preds[i] = predictMurine(dft[i, :]; kwargs...)
    end
    dft."Predict" = preds
    if any(dft."Predict" .> 0.0)
        dft."Predict" ./= geomean(dft."Predict"[dft."Predict" .> 0.0]) / geomean(dft."Value"[dft."Predict" .> 0.0])
    end

    return dft
end

function predictMurine(
    c::Chains = runMurineMCMC(),
    df::AbstractDataFrame = importMurineInVitro();
    Kav = murineKavDist(; regularKav = true),
    kwargs...,
)
    Rtot = [median(c["Rtot[$i]"].data) for i = 1:4]
    Rtotd = Dict([murineFcgR[ii] => Rtot[ii] for ii = 1:length(Rtot)])
    f33 = median(c["f33"].data)
    KxStar = median(c["KxStar"].data)
    return predictMurine(df; Kav = Kav, KxStar = KxStar, recepExp = Rtotd, f33 = f33, kwargs...)
end

@model function murineFit(df, values; Kavd::Union{Nothing, AbstractDataFrame} = nothing)
    # Rtot sample
    Rtot = Vector(undef, length(murineFcgR))
    for ii in eachindex(Rtot)
        ref = InVitroMurineRcpExp[murineFcgR[ii]]
        Rtot[ii] ~ inferLogNormal(ref, ref * 1e2)
        # Treat receptor amount as unknown with a wide prior
    end
    Rtotd = Dict([murineFcgR[ii] => Rtot[ii] for ii = 1:length(Rtot)])

    if Kavd === nothing
        # fit Kav
        Kav_dist = Matrix(murineKavDist()[:, Not("IgG")])
        Kavd = murineKavDist(; regularKav = true)
        Kavd = Kavd[Kavd."IgG" .!= "IgG3", :]
        Kav = Matrix(undef, size(Kav_dist)...)
        for ii in eachindex(Kav)
            Kav[ii] ~ Kav_dist[ii]
        end
        Kavd[!, Not("IgG")] = Kav
    end

    f33 ~ f33Dist
    KxStar ~ KxStarDist

    # fit predictions
    if all(Rtot .>= 0.0) && all(0.0 .<= Matrix(Kavd[!, Not("IgG")]) .< Inf)
        df = predictMurine(deepcopy(df); recepExp = Rtotd, Kav = Kavd, KxStar = KxStar, f33 = f33)
    else
        df = deepcopy(df)
        df."Predict" .= Inf
    end

    stdv = std(log.(df."Predict") - log.(values))
    values ~ MvLogNormal(log.(df."Predict"), stdv * I)
    nothing
end

function runMurineMCMC(fname = "murineNUTSdepfit_0505.dat"; mcmc_iter = 1_000)
    if isfile(fname)
        return deserialize(fname)
    end
    df = importMurineInVitro()
    m = murineFit(df, df."Value")
    # use MAP estimate as starting point
    opts = Optim.Options(iterations = 200, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    c = sample(m, NUTS(), mcmc_iter, init_params = opt.values.array)
    f = serialize(fname, c)
    return c
end


function plot_murineMCMC_predict(c::Chains = runMurineMCMC(), df = importMurineInVitro(); title = nothing, kwargs...)
    p = extractMCMC(c; murine = true)
    title = (title === nothing) ? "Murine fitting with NUTS (median)" : title

    ndf = predictMurine(df; recepExp = p["Rtot"], Kav = p["Kav"], KxStar = p["KxStar"], f33 = p["f33"])
    return plotPredvsMeasured(ndf; xx = "Value", yy = "Predict", color = "Receptor", shape = "Subclass", R2pos = (0, -1), title = title, kwargs...)
end

function MAPmurineLikelihood(df = importMurineInVitro())
    m = murineFit(df, df."Value")
    opts = Optim.Options(iterations = 200, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    x = opt.values

    Rtot = Dict([rcp => x[Symbol("Rtot[$i]")] for (i, rcp) in enumerate(murineFcgR)])
    Kav = murineKavDist(; regularKav = true, retdf = true)
    Kav[!, Not("IgG")] = reshape([x[Symbol("Kav[$i]")] for i = 1:12], 3, 4)

    ndf = predictMurine(df; Kav = Kav, recepExp = Rtot)
    pl = plotPredvsMeasured(
        ndf;
        xx = "Value",
        yy = "Predict",
        color = "Receptor",
        shape = "Subclass",
        R2pos = (-0.5, -1.2),
        title = "Murine MAP fitting results",
    )
    return pl, [Rtot, Kav]
end

function validateMurineInVitro(c::Chains = fitLeukocyteMCMC(); mcmc_iter = 1_000)
    df = importMurineInVitro()
    opts = Optim.Options(iterations = 200, show_every = 10, show_trace = true)

    Kav_old = importKav(; murine = true, retdf = true)
    m_old = murineFit(df, df."Value"; Kavd = Kav_old)
    opt_old = optimize(m_old, MAP(), LBFGS(; m = 20), opts)
    c_old = sample(m_old, NUTS(), mcmc_iter, init_params = opt_old.values.array)

    Kav_new = extractMCMC(c; murine = true)["Kav"]
    m_new = murineFit(df, df."Value"; Kavd = Kav_new)
    opt_new = optimize(m_new, MAP(), LBFGS(; m = 20), opts)
    c_new = sample(m_new, NUTS(), mcmc_iter, init_params = opt_new.values.array)

    pl1 = plotPredvsMeasured(
        predictMurine(c_old, df; Kav = Kav_old);
        xx = "Value",
        yy = "Predict",
        color = "Receptor",
        shape = "Subclass",
        R2pos = (-1.5, 0.0),
        title = "Murine in vitro binding prediction\nwith documented affinities",
    )
    pl2 = plotPredvsMeasured(
        predictMurine(c_new, df; Kav = Kav_new);
        xx = "Value",
        yy = "Predict",
        color = "Receptor",
        shape = "Subclass",
        R2pos = (-1.5, 0.2),
        title = "Murine in vitro binding prediction\nwith updated affinities",
    )

    #pp = plotGrid((1, 2), [pl1, pl2])
    #draw(PDF("figure3CHO.pdf", 7inch, 3inch), pp)
    return pl1, pl2
end
