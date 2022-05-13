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

const InVitroMurineRcpExp = Dict("FcgRI" => 1e6, "FcgRIIB" => 1e6, "FcgRIII" => 1e6, "FcgRIV" => 1e6)

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
        Kav_dist = Matrix(importKavDist(; murine = true, regularKav = false)[:, Not("IgG")])
        Kavd = importKavDist(; murine = true, regularKav = true)
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
        df = predMix(deepcopy(df); Rtot = Rtotd, Kav = Kavd, KxStar = KxStar, fs = [4, f33])
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


function MAPmurineLikelihood(df = importMurineInVitro())
    m = murineFit(df, df."Value")
    opts = Optim.Options(iterations = 200, show_every = 10, show_trace = true)
    opt = optimize(m, MAP(), LBFGS(; m = 20), opts)
    x = opt.values

    Rtot = Dict([rcp => x[Symbol("Rtot[$i]")] for (i, rcp) in enumerate(murineFcgR)])
    Kav = importKavDist(; murine = true, regularKav = true, retdf = true)
    Kav[!, Not("IgG")] = reshape([x[Symbol("Kav[$i]")] for i = 1:12], 3, 4)

    ndf = predMix(df; Kav = Kav, Rtot = Rtot)
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

function validateMurineInVitro(c::Chains = rungMCMC("leukNUTSfit_0509.dat"); mcmc_iter = 1_000)
    df = importMurineInVitro()
    opts = Optim.Options(iterations = 200, show_every = 10, show_trace = true)

    Kav_old = importKav(; murine = true, retdf = true)
    m_old = gmodel(df, df."Value"; dat = :mCHO, Kavd = Kav_old)
    opt_old = optimize(m_old, MAP(), LBFGS(; m = 20), opts)
    c_old = sample(m_old, NUTS(), mcmc_iter, init_params = opt_old.values.array)

    Kav_new = extractMCMC(c; dat = :mLeuk)["Kav"]
    m_new = gmodel(df, df."Value"; dat = :mCHO, Kavd = Kav_new)
    opt_new = optimize(m_new, MAP(), LBFGS(; m = 20), opts)
    c_new = sample(m_new, NUTS(), mcmc_iter, init_params = opt_new.values.array)

    pl1 = plotMCMCPredict(c_old, df; dat = :mCHO, Kav = Kav_old, 
        R2pos = (-1.5, 0.0), title = "Murine in vitro binding prediction\nwith documented affinities")
    pl2 = plotMCMCPredict(c_new, df; dat = :mCHO, Kav = Kav_new, 
        R2pos = (-1.5, 0.2), title = "Murine in vitro binding prediction\nwith updated affinities")
    return pl1, pl2
end
