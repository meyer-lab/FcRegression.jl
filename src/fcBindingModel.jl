function polyfc_ActV(L0, KxStar, f, Rtot::Array, IgGC::Array, Kav::AbstractMatrix, Rbound::Bool = false; Mix = true)
    """
    Input:
    Rtot: nr * nct matrix, nct = # cell types
    IgGC: ni * nset matrix, nset = # IgGC combinations
    Mix: Specifies that output will be the response to a single IgG input if mix = false
    Output:
    Matrix of size nct * nset filled with Rmulti_n/Rbound (depends on if Rbound is provided)
    """
    ansType = promote_type(typeof(L0), typeof(KxStar), typeof(f), eltype(Rtot), eltype(IgGC))
    if !Rbound
        res = Array{ansType}(undef, size(Rtot, 2), size(Rtot, 1), size(IgGC, 2))
    else
        @assert ndims(Rtot) == 1
        res = Vector{ansType}(undef, size(IgGC, 2))
    end
    @assert size(Rtot, 1) == size(Kav, 2)
    @assert size(IgGC, 1) == size(Kav, 1)

    if Mix
        L0 = ones(size(IgGC, 2)) * L0
    else
        L0 = (range(0, stop = 1, length = size(IgGC, 2))) * L0
    end
    for iset = 1:size(IgGC, 2)
        if Rbound
            res[iset] = polyfc(L0[iset], KxStar, f, Rtot, IgGC[:, iset], Kav).Rbound
        else
            for ict = 1:size(Rtot, 2)
                res[ict, :, iset] .= polyfc(L0[iset], KxStar, f, Rtot[:, ict], IgGC[:, iset], Kav).Rmulti_n
            end
        end
    end
    return res
end
