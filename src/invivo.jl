function invivo_mix_predict(dfr::DataFrameRow; L0 = 1e-9, f = 4, KxStar = KxConst, murine = true)
    IgGs = murine ? murineIgG : humanIgG
    IgGC = zeros(size(IgGs))
    IgGXname = (dfr."subclass_1" == "IgG2c") ? "IgG2a" : dfr."subclass_1"
    IgGYname = (dfr."subclass_2" == "IgG2c") ? "IgG2a" : dfr."subclass_2"
    IgGC[IgGXname .== IgGs] .= dfr."%_1"
    IgGC[IgGYname .== IgGs] .= dfr."%_2"
    if sum(IgGC) <= 0.0
        L0 = 0.0
        IgGC[1] = 1
    else
        IgGC ./= sum(IgGC)
    end

    Kav = importKav(; murine = murine)
    Rtot = importRtot(; murine = murine)
    ActI = murine ? murineActI : humanActI
    return polyfc_ActV(L0, KxStar, f, Rtot, IgGC, Kav, false)[:, :, 1] * ActI
end

function invivo_prediction(df::DataFrame; L0 = 1e-9, f = 4, KxStar = KxConst)
    fcPredict = cat([invivo_mix_predict(dfr; L0 = L0, f = f, KxStar = KxConst) for dfr in eachrow(df)]...; dims = 2)'
    res = regResult("ITP"; L0 = 1e-8, f = 10, murine = true)[1]
    Ypred = regPred(Matrix(fcPredict), res; exp_method = true)
    df[!, "Predicted"] = Ypred
    return df
end

""" Make a descriptive name for an IgG mixture dose"""
function generate_mixture_name(dfr::DataFrameRow)
    str = ""
    if dfr."%_1" > 0
        str *= dfr."subclass_1" * " " * string(dfr."%_1") * "%"
    end
    if dfr."%_2" > 0
        if str > ""
            str *= ", "
        end
        str *= dfr."subclass_2" * " " * string(dfr."%_2") * "%"
    end
    if str <= ""
        str *= "PBS"
    end
    return str
end
