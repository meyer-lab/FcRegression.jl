
function polyfc_ActV(L0, KxStar, f, Rtot::Array, IgGC::Array, Kav::AbstractMatrix, ActI::Vector; Mix = true)
    """
    Input:
    Rtot: nr * nct matrix, nct = # cell types
    IgGC: ni * nset matrix, nset = # IgGC combinations
    Mix: Specifies that output will be the response to a single IgG input if mix = false 
    Output:
    Matrix of size nct * nset filled with ActV
    """
    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))
    res = Matrix{ansType}(undef, size(Rtot, 2), size(IgGC, 2))
    if Mix
        L0 = ones(size(res, 2)) * L0
    else
        L0 = (range(0, stop = 1, length = size(res, 2))) * L0
    end
    for ict = 1:size(res, 1)
        for iset = 1:size(res, 2)
            res[ict, iset] = polyfc(L0[iset], KxStar, f, Rtot[:, ict], IgGC[:, iset], Kav, ActI).ActV
        end
    end
    return res
end


function polyfcm_ActV(IgG::AbstractMatrix, KxStar, f, Rtot, Kav, ActI)
    """
    Input:
    Rtot: nr * nct matrix, nct = # cell types
    IgG: ni * nset matrix, nset = # IgG combinations
    Output:
    Matrix of size nct * nset filled with ActV
    """
    ansType = promote_type(typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgG))
    res = Matrix{ansType}(undef, size(Rtot, 2), size(IgG, 2))
    for ict = 1:size(res, 1)
        for iset = 1:size(res, 2)
            res[ict, iset] = polyfcm(KxStar, f, Rtot[:, ict], IgG[:, iset], Kav, ActI).ActV
        end
    end
    return res
end
