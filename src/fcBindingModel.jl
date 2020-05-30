
function polyfc_ActV(L0, KxStar, f, Rtot::Array, IgGC::Array, Kav::AbstractMatrix, ActI::Vector; Mix = true)
    """
    Input:
    Rtot: nr * nct matrix, nct = # cell types
    IgGC: ni * nset matrix, nset = # IgGC combinations
    Output:
    Matrix of size nct * nset filled with ActV
    """
    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))
    res = Matrix{ansType}(undef, size(Rtot, 2), size(IgGC, 2))
    if Mix
        for ict = 1:size(res, 1)
            for iset = 1:size(res, 2)
                res[ict, iset] = polyfc(L0, KxStar, f, Rtot[:, ict], IgGC[:, iset], Kav, ActI).ActV
            end
        end
    else
        L = (range(0, stop = 1, length = size(res, 2)))*L0
        for ict = 1:size(res, 1)
            for iset = 1:size(res, 2)
                res[ict, iset] = polyfc(L[iset], KxStar, f, Rtot[:, ict], IgGC[:, iset], Kav, ActI).ActV
            end
        end
    end
    return res
end
