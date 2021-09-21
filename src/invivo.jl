function invivo_mix_predict(dfr::DataFrameRow; L0 = 1e-9, f = 4, KxStar = KxConst)
    IgGC = zeros(size(murineIgG))
    IgGXname = (dfr."subclass_1" == "IgG2c") ? "IgG2a" : dfr."subclass_1"
    IgGYname = (dfr."subclass_2" == "IgG2c") ? "IgG2a" : dfr."subclass_2"
    IgGC[IgGXname .== murineIgG] .= dfr."%_1"
    IgGC[IgGYname .== murineIgG] .= dfr."%_2"
    if sum(IgGC) <= 0.0
        L0 = 0.0
        IgGC[1] = 1.0
    else
        IgGC ./= sum(IgGC)
    end

    Kav = importKav(; murine = true)
    Rtot = importRtot(; murine = true)
    return polyfc_ActV(L0, KxStar, f, Rtot, IgGC, Kav, false)
end

function invivo_prediction(df::DataFrame; L0 = 1e-9, f = 4, KxStar = KxConst)
    fcPredict = cat([invivo_mix_predict(dfr; L0 = L0, f = f, KxStar = KxConst) for dfr in eachrow(df)]...; dims = 3)
    res, _, _, _ = regressionResult("ITP"; L0 = 1e-9, f = 10, murine = true)
    Xmat, Ypred = regressionPred(fcPredict, nothing, res.cellWs, res.ActI; showXmat = true, murine = true)
    df[!, "Predicted"] = exponential(Ypred)
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
