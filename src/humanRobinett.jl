function figure2d()
    df = averageData(loadMixData("robinett/Luxetal2013-Fig2BRef.csv"))
    df = predictMix(df; recepExp = RobMeasuredRecepExp)
    df."Adjusted" = df."Value" .* (geocmean(df."Predict") / geocmean(df."Value"))
    draw(SVG("figure2d.svg", 1300px, 600px), plotGrid((1, 2), [nothing, plotPredvsMeasured(df)]))
end

function importRobinett()
    df = CSV.File(joinpath(dataDir, "robinett/Luxetal2013-Fig2Bmod.csv"), delim = ",", comment = "#") |> DataFrame
    for i = 1:4
        cn = "Replicate $i"
        df[!, cn] ./= geomean(df[Not(ismissing.(df[!, cn])), cn])
    end
    df = dropmissing(stack(df, Not(["Cell", "Antibody", "Valency"])))
    rename!(df, ["variable" => "Experiment", "value" => "Value"])
    rename!(df, ["Antibody" => "subclass_1"])
    df[!, "%_1"] .= 1.0
    df[!, "subclass_2"] .= "None"
    df[!, "%_2"] .= 0.0
    return sort!(df, ["Valency", "Cell", "subclass_1", "subclass_2", "Experiment", "%_2"])
end

""" Fit everything but affinities for Robinett data with MCMC, compare new and old affinities"""
function validateRobinett(fname = "MCMC_robinett_0505.dat", c = runMCMC(); mcmc_iter = 1_000)
    local c_old, c_new
    if isfile(fname)
        c_old, c_new = deserialize(fname)
        df = importRobinett()
    else
        df = importRobinett()
        opts = Optim.Options(iterations = 500, show_every = 10, show_trace = true)

        Kav_old = importKav(; murine = false, invitro = true, retdf = true)
        m_old = sfit(df, df."Value"; robinett = true, Kavd = Kav_old)
        opt_old = optimize(m_old, MAP(), LBFGS(; m = 20), opts)
        c_old = sample(m_old, NUTS(), mcmc_iter, init_params = opt_old.values.array)

        _, Kav_new, _ = extractMCMCresults(c)
        m_new = sfit(df, df."Value"; robinett = true, Kavd = Kav_new)
        opt_new = optimize(m_new, MAP(), LBFGS(; m = 20), opts)
        c_new = sample(m_new, NUTS(), mcmc_iter, init_params = opt_new.values.array)

        f = serialize(fname, [c_old, c_new])
    end

    pl1 = MCMC_params_predict_plot(c_old, df; xx = "Value", yy = "Predict", title = "Robinett with documented affinities", R2pos = (0, -1.5))
    pl2 = MCMC_params_predict_plot(c_new, df; xx = "Value", yy = "Predict", title = "Robinett with updated affinities", R2pos = (0, -1.5))

    pp = plotGrid((1, 2), [pl1, pl2])
    draw(PDF("figure2robinett.pdf", 7inch, 3inch), pp)
    return pl1, pl2
end
