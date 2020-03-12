using Optim
import MLBase.LOOCV
import Statistics: mean, std
import Distributions: cdf, Exponential
import StatsBase: sample, mode
using GLM
using NNLS

exponential(X::Matrix, p::Vector) = cdf.(Exponential(), X * p)
inv_exponential(y::Real) = -log(1 - y)

function regGenData(df; L0, f, murine = true, retdf = false, component = cellTypes)
    df = copy(df)
    Rtot = importRtot(; murine = murine)
    FcRecep = murine ? murineFcgR : humanFcgR
    ActI = murine ? murineActI : humanActI

    if :Concentration in names(df)
        df[!, :Concentration] .*= L0
    else
        insertcols!(df, 3, :Concentration => L0)
    end

    X = DataFrame(repeat([Float64], length(cellTypes)), cellTypes)
    for i = 1:size(df, 1)
        Kav = convert(Vector{Float64}, df[i, FcRecep])
        Kav = reshape(Kav, 1, :)
        push!(X, polyfc_ActV(df[i, :Concentration], KxConst, f, Rtot, [1.0], Kav, ActI))
    end

    # any extra column (C1q or Neutralization) will be taken directly as an extra term in regression
    for extra_col in setdiff(names(df), [:Condition; :Background; :Concentration; :Target; FcRecep])
        X[!, extra_col] = df[!, extra_col]
    end

    X[df[:, :Background] .== "NeuKO", :Neu] .= 0.0
    X[df[:, :Background] .== "ncMOKO", :ncMO] .= 0.0
    Y = df[!, :Target]

    @assert all(isfinite.(Matrix(X)))
    @assert all(isfinite.(Y))

    if retdf
        return (X, Y)
    else
        return (Matrix{Float64}(X), Y)
    end
end


function quadratic_loss(Y0::Vector, Y::Vector)
    return Distances.sqeuclidean(Y0, Y)
end

function proportion_loss(p::Vector, Y::Vector)
    res = sum((Y .- p) .^ 2 ./ (p .* (1 .- p) .+ 0.25))
    @assert all(isfinite.(res))
    return res
end

function loss_wL0f(df, ps::Vector{T}, lossFunction::Function)::T where {T <: Real}
    (X, Y) = regGenData(df; L0 = 10.0^ps[1], f = ps[2])
    Y0 = exponential(X, ps[3:end])
    return lossFunction(Y0, Y)
end

function fitRegression(df, lossFunction::Function = proportion_loss; wL0f = false)
    ## this method only supports expoential distribution due to param choice

    (X, Y) = regGenData(df; L0 = 1.0e-9, f = 4, retdf = false)
    Np = size(X, 2)
    if wL0f
        fitMethod = (ps) -> loss_wL0f(df, ps, lossFunction)
    else
        fitMethod = (ps) -> lossFunction(exponential(X, ps), Y)
    end
    g! = (G, ps) -> ForwardDiff.gradient!(G, fitMethod, ps)

    p_init = nnls(X, inv_exponential.(Y))
    p_lower = zeros(Float64, Np)
    p_upper = maximum(p_init) .* ones(Float64, Np)

    if wL0f
        p_init = vcat(-9, 4, p_init)
        p_lower = vcat(-16, 1, p_lower)
        p_upper = vcat(-7, 6, p_upper)
    end

    fit = optimize(fitMethod, g!, p_lower, p_upper, p_init, Fminbox())
    if !Optim.converged(fit)
        @warn "Fitting did not converge"
    end
    return p_init
end

function LOOCrossVal(dataType, lossFunction::Function; wL0f = false)
    df = importDepletion(dataType)
    n = size(df, 1)
    fitResults = Vector(undef, n)
    LOOindex = LOOCV(n)
    for (i, idx) in enumerate(LOOindex)
        fitResults[i] = fitRegression(df[idx, :], lossFunction; wL0f = wL0f)
    end
    return fitResults
end

function bootstrap(dataType, lossFunction::Function; nsample = 100, wL0f = false)
    df = importDepletion(dataType)
    n = size(df, 1)
    fitResults = Vector(undef, nsample)
    for i = 1:nsample
        for j = 1:5
            fit = try
                fitRegression(df[sample(1:n, n, replace = true), :], lossFunction; wL0f = wL0f)
            catch e
                @warn "This bootstrapping set failed at fitRegression"
                nothing
            end
            if fit != nothing
                fitResults[i] = fit
                break
            end
        end
    end
    @assert all(isa.(fitResults, Array{<:Real}))
    return fitResults
end


function CVResults(dataType, lossFunction::Function = proportion_loss; L0 = 1e-9, f = 4)
    df = importDepletion(dataType)
    fit_w = fitRegression(df, lossFunction)
    loo_out = LOOCrossVal(dataType, lossFunction)
    btp_out = bootstrap(dataType, lossFunction)

    (X, Y) = regGenData(df; L0 = L0, f = f, retdf = true)
    components = names(X)
    @assert length(fit_w) == length(components)

    odf = df[!, [:Condition, :Background]]
    odf[!, :Concentration] .= (:Concentration in names(df)) ? (df[!, :Concentration] .* L0) : L0
    odf[!, :Y] = Y
    odf[!, :Fitted] = exponential(Matrix(X), fit_w)
    odf[!, :LOOPredict] = vcat([exponential(Matrix(X[[i], :]), loo_out[i]) for i = 1:length(loo_out)]...)

    selected = (odf[!, :Background] .== "wt") .& (odf[!, :Concentration] .== mode(odf[!, :Concentration]))
    X = Matrix(X[selected, :])

    effects = X .* fit_w'
    btp_ws = cat([X .* a' for a in btp_out]..., dims = (3))
    btp_std = dropdims(std(btp_ws; dims = 3), dims = 3)
    @assert size(effects) == size(btp_std)

    wdf = DataFrame(Weight = vec(effects), BtpStdev = vec(btp_std))
    wdf[!, :Condition] = vec(repeat(odf[selected, :Condition], 1, length(fit_w)))
    wdf[!, :Component] = vec(repeat(reshape(components, 1, :), size(effects, 1), 1))
    @assert size(wdf, 1) == size(X, 1) * length(fit_w)

    return (fit_w, odf, wdf)
end



function regGenHumanized(df; L0 = 1e-9, f = 4, murine = false)
    X = []
    for dfr in eachrow(df)
        genotype = :Genotype in names(dfr) ? dfr.:Genotype : "HIV"
        Lconc = :Concentration in names(dfr) ? dfr.:Concentration * L0 : L0

        IgGnames = murine ? murineIgG : humanIgG
        IgGC = zeros(length(IgGnames))
        IgGC[IgGnames .== dfr.:Condition] .= 1
        @assert sum(IgGC) == 1

        Rtot = importRtot(; murine = murine, genotype = genotype)
        Kav = importKav(; murine = murine, c1q = true, retdf = true)
        FcRecep = murine ? murineFcgR : humanFcgR
        ActI = murine ? murineActI : humanActI

        res = polyfc_ActV(Lconc, KxConst, f, Rtot, IgGC, Matrix(Kav[!, FcRecep]), ActI)
        effects = Dict(cellTypes .=> res)
        effects[:C1q] = Kav[!, :C1q]' * IgGC .* Lconc
        push!(X, effects)
    end

    X = DataFrame(X)
    Y = df[!, :Target]

    return (X, Y)
end
