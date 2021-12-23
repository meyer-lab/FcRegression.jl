""" Fitting mixture measurements """

using Optim

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
    Kav::DataFrame = importKav(; murine = false, retdf = true),
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
    Kav::DataFrame = importKav(; murine = false, retdf = true),
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
        return deepcopy(df[!, ["IgG"; humanFcgR]])
    else
        return deepcopy(Matrix{Float64}(df[!, humanFcgR]))
    end
end
