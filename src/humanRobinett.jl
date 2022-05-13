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

""" Fit everything but affinities for Robinett data with MCMC, compare new and old affinities"""
function validateRobinett(fname = "MCMC_robinett_0505.dat", c = rungMCMC("humanNUTSfit_0505.dat"); mcmc_iter = 1_000)
    local c_old, c_new
    df = importRobinett()
    Kav_old = importKavDist(; murine = false, regularKav = true, retdf = true)
    Kav_new = extractMCMC(c; dat = :hCHO)["Kav"]

    if isfile(fname)
        c_old, c_new = deserialize(fname)
    else
        opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)

        m_old = gmodel(df, df."Value"; dat = :hRob, Kavd = Kav_old)
        opt_old = optimize(m_old, MAP(), LBFGS(; m = 20), opts)
        c_old = sample(m_old, NUTS(), mcmc_iter, init_params = opt_old.values.array)

        m_new = gmodel(df, df."Value"; dat = :hRob, Kavd = Kav_new)
        opt_new = optimize(m_new, MAP(), LBFGS(; m = 20), opts)
        c_new = sample(m_new, NUTS(), mcmc_iter, init_params = opt_new.values.array)

        f = serialize(fname, [c_old, c_new])
    end

    pl1 = plotMCMCPredict(c_old, df; dat = :hRob, Kav = Kav_old,
        title = "Robinett with documented affinities", R2pos = (-0.5, -2))
    pl2 = plotMCMCPredict(c_new, df; dat = :hRob, Kav = Kav_new,
        title = "Robinett with updated affinities", R2pos = (-0.5, -2))
    
    return pl1, pl2
end
