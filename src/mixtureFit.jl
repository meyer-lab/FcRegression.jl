""" Fitting mixture measurements """

using Optim

function mixturePredictions(df = loadMixData();
    Rtot = measuredRecepExp,
    Kav = importKav(; murine = false, invitro = true, retdf = true),
    KxStar = KxConst,
    vals = [4.0, 33.0],
)
    try
        df[!, "NewValency"] .= vals[1]
        df[df."Valency" .== 33, "NewValency"] .= vals[2]
    catch e
        println("Rtot, ", eltype(Rtot))
        println("KxStar, ", KxStar)
        println("vals, ", vals)
        println("df.NewVal, ", eltype(df."NewValency"))
        rethrow(e)
    end
    #println("** in mixturePredictions()")
    #println(Rtot, Kav, KxStar, vals)
    ndf = predictMix(df; recepExp = Rtot, KxStar = KxStar, Kav = Kav)
    #valconv1 = sum(log.(ndf[ndf."Valency" .== 4, "Value"])) / sum(log.(ndf[ndf."Valency" .== 4, "Predict"]))
    #valconv2 = sum(log.(ndf[ndf."Valency" .== 33, "Value"])) / sum(log.(ndf[ndf."Valency" .== 33, "Predict"]))
    #ndf[ndf."Valency" .== 4, "Predict"] .*= valconv1
    #ndf[ndf."Valency" .== 33, "Predict"] .*= valconv2
    return ndf
end

function fit_goodness(ndf; deviation = 0.01)
    # log likelihood of prediction vs measurements
    # deviation < 1.0
    llike = 0
    for ii = 1:size(ndf, 1)
        lmval = ndf[ii, "Value"]
        llike += log(pdf(Normal(lmval, lmval * deviation), log(ndf[ii, "Predict"])))
    end
    return llike
end

function Rtot_prior(Rtot::Dict)
    # return log f(R1)f(R2)...f(Rn), input regular Ri w/o log
    # to speed up, order (lexicographical) is not checked here, but it is important!!
    dists = importInVitroRtotDist()
    return sum([log.(pdf(dists[ii], log(Rtot[humanFcgRiv[ii]]))) for ii = 1:length(Rtot)])
end

function Kav_prior(Kpropose::DataFrame)
    # return log f(Ka1) f(Ka2) ... f(Kan), input regular Kai w/o log
    Kav = importKavDist(; inflation = 0.1)
    @assert size(Kav) == size(Kpropose)
    Kav = reshape(log.(Matrix(Kav[:, Not("IgG")])), :)
    Kpropose = reshape(log.(Matrix(Kpropose[:, Not("IgG")])), :)
    Kpdf = (x_dist, x_test) -> pdf(x_dist, x_test)
    return sum(log.([Kpdf(Kav[ii], Kpropose[ii]) for ii = 1:length(Kav)]))
end

function valency_prior(vals; vstd = 0.2)
    # return log f(Val4) f(Val33), input regular val w/o log
    return log(pdf(Normal(log(4), vstd), log(vals[1]))) + log(pdf(Normal(log(33), vstd), log(vals[2])))
end

KxStar_prior = x -> log(pdf(Normal(log(KxConst), 4.37), log(x)))    # eyeballed from Robinett Fig. 2d

function dismantle_x0(x)
    # order: Rtot, vals, KxStar, Kav
    x = deepcopy(x)
    @assert all(x[1:length(humanFcgRiv)] .> 0.0) "In dismantle(): Rtot contains negative"
    Rtot = Dict([humanFcgRiv[ii] => popfirst!(x) for ii = 1:length(humanFcgRiv)])
    vals = [popfirst!(x), popfirst!(x)]
    @assert all(vals .> 0.0) "In dismantle(): vals contains negative"
    KxStar = popfirst!(x)
    @assert KxStar .> 0.0 "In dismantle(): KxStar is negative"
    @assert length(x) == length(humanIgG) * length(humanFcgRiv)
    @assert all(x .> 0.0) "In dismantle(): Kav contains negative"
    Kav = deepcopy(importKav(; murine = false, invitro = true, retdf = true))
    Kav[:, Not("IgG")] .= 0.0
    Kav[!, Not("IgG")] = convert.(eltype(x), Kav[!, Not("IgG")])
    Kav[:, Not("IgG")] = reshape(x, length(humanIgG), length(humanFcgRiv))
    return Rtot, vals, KxStar, Kav
end

function assemble_x0(
    Rtot::Dict = measuredRecepExp, 
    vals::Vector = [4.0, 33.0], 
    KxStar = KxConst, 
    Kav::DataFrame = importKav(; murine = false, invitro = true, retdf = true)
)
    x = [Rtot[rr] for rr in humanFcgRiv]
    push!(x, vals...)
    push!(x, KxStar)
    push!(x, reshape(Matrix(Kav[:, Not("IgG")]), :)...)
    return x
end

function totalLikelihood(x, df = loadMixData(; discard_small = true); deviation = 0.1)
    Rtot, vals, KxStar, Kav = dismantle_x0(x)
    ndf = mixturePredictions(df; Rtot = Rtot, Kav = Kav, KxStar = KxStar, vals = vals)
    lik = Rtot_prior(Rtot) + Kav_prior(Kav) + valency_prior(vals) + KxStar_prior(KxStar)
    lik += fit_goodness(ndf; deviation = deviation)
    return lik
end


function MLELikelihood()
    x0 = log.(FcRegression.assemble_x0());
    f = lx -> -FcRegression.totalLikelihood(exp.(lx); deviation = 0.01)       # minimize the negative of likelihood
    ## TODO: write Optim function

    # GD: 7.243698e+04 (reach maxiter)
    # BFGS: 7.343904e+04
    # LBFGS: 7.201175e+04
    opt = optimize(f, x0, LBFGS(), Optim.Options(iterations = 100, show_trace = true); autodiff = :forward)

    Rtot, vals, KxStar, Kav = FcRegression.dismantle_x0(exp.(opt.minimizer))
    ndf = FcRegression.mixturePredictions(FcRegression.averageMixData(FcRegression.loadMixData()); Rtot = Rtot, Kav = Kav, KxStar = KxStar, vals = vals);
    FcRegression.plotPredvsMeasured(ndf; xx = "Value")
end





















#=

function mixSqLoss(df; logscale = true)
    # square root differences of model prediction and adjusted measurements
    if "Adjusted" in names(df)
        adj = logscale ? log.(df."Adjusted") : df."Adjusted"
    else
        adj = logscale ? log.(df."Value") : df."Value"
    end
    pred = logscale ? log.(df."Predict") : df."Predict"
    return sum((adj .- pred) .^ 2)
end


function fitMixFunc(
    x::Vector,
    df;
    fitRVX = true,
    recepExp::Union{Vector, Dict} = measuredRecepExp,
    vals = [4, 33],
    KxStar = KxConst,
    fitKav = false,
    Kav::DataFrame = importKav(; murine = false, invitro = true, retdf = true),
)
    ## assemble numbers to a vector for optimization
    # order: log(Rtot), log(valency), log(Kx*), log(Kav)

    @assert fitRVX | fitKav
    cells = sort(unique(df."Cell"))
    @assert length(x) == (fitRVX ? length(cells) + 3 : 0) + (fitKav ? size(Kav)[1] * (size(Kav)[2] - 1) : 0)
    x = exp.(x)

    if fitRVX
        recepExp = Dict([cell => x[ii] for (ii, cell) in enumerate(cells)])
        # adjust the valency, given there are only 4 and 33 originally
        df[!, "NewValency"] .= x[(length(cells) + 1)]
        df[df."Valency" .== 33, "NewValency"] .= x[(length(cells) + 2)]
        KxStar = x[(length(cells) + 3)]
    else
        df[!, "NewValency"] .= vals[1]
        df[df."Valency" .== 33, "NewValency"] .= vals[2]
        if recepExp isa Vector
            recepExp = Dict([cell => recepExp[ii] for (ii, cell) in enumerate(cells)])
        end
    end
    if fitKav
        Kav[!, Not("IgG")] .= 0.0
        Kav[!, Not("IgG")] = reshape((fitRVX ? x[(length(cells) + 4):end] : x), (size(Kav)[1], :))
    end

    return predictMix(df; recepExp = recepExp, KxStar = KxStar, Kav = Kav)
end

function fitMixMaster(
    df = loadMixData();
    fitRVX = true,
    recepExp = measuredRecepExp,
    vals = [4, 33],
    KxStar = KxConst,
    fitKav = false,
    Kav::DataFrame = importKav(; murine = false, invitro = true, retdf = true),
    show_trace = false,
)
    # order: log(Rtot), log(valency), log(Kx*), log(Kav)
    x0 = Vector{Float64}()

    # x0 for Rtot, valency, Kx*
    cells = sort(unique(df."Cell"))
    if fitRVX
        @assert recepExp isa Dict
        x0 = [log(recepExp[cell]) for cell in cells]
        append!(x0, log.(vals))
        append!(x0, log(KxConst))
    end

    # x0 for Kav
    if fitKav
        ### TODO: add options to only fit some Kav
        x0aff = reshape(Matrix(Kav[!, Not("IgG")]), (:))
        x0aff[x0aff .<= 0.0] .= 1.0
        append!(x0, log.(x0aff))
    end

    # set bounds for optimization
    x_ub = x0 .+ 3.0
    x_lb = x0 .- 3.0
    if fitRVX
        x_ub[(length(cells) + 1):(length(cells) + 2)] .= x0[(length(cells) + 1):(length(cells) + 2)] .+ 0.01   # valency ub is the original valency
        x_lb[(length(cells) + 1):(length(cells) + 2)] .= x0[(length(cells) + 1):(length(cells) + 2)] .- 1.0    # valency lb is exp(1) below
    end
    println(x0, x_ub, x_lb)
    f = x -> mixSqLoss(fitMixFunc(x, df; fitRVX = fitRVX, recepExp = recepExp, vals = vals, KxStar = KxStar, fitKav = fitKav, Kav = Kav))

    dfc = TwiceDifferentiableConstraints(x_lb, x_ub)
    res = optimize(f, dfc, x0, IPNewton(), Optim.Options(iterations = 100, show_trace = show_trace); autodiff = :forward)
    ndf =
        fitMixFunc(res.minimizer, averageMixData(df); fitRVX = fitRVX, recepExp = recepExp, vals = vals, KxStar = KxStar, fitKav = fitKav, Kav = Kav)
    if fitKav
        Kav[!, Not("IgG")] .= 0.0
        Kav[!, Not("IgG")] = reshape(exp.((fitRVX ? res.minimizer[(length(cells) + 4):end] : res.minimizer)), (size(Kav)[1], :))
        nKav = unstack(stack(Kav, Not("IgG")), "variable", "IgG", "value")
        CSV.write(joinpath(dataDir, "fitted_human_new_affinity.csv"), nKav)
        return res, ndf, Kav
    end
    return res, ndf
end

function loadFittedKav(; retdf = true)
    df = CSV.File(joinpath(dataDir, "fitted_human_new_affinity.csv"), comment = "#") |> DataFrame
    df = stack(df; variable_name = "IgG", value_name = "Kav")
    df = unstack(df, "variable", "Kav")
    df = df[in(humanIgG).(df.IgG), :]
    if retdf
        return deepcopy(df[!, ["IgG"; humanFcgRiv]])
    else
        return deepcopy(Matrix{Float64}(df[!, humanFcgRiv]))
    end
end

=#