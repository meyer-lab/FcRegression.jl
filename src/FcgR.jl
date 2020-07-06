module FcgR
using LinearAlgebra
using ForwardDiff
import Distances
using Gadfly
using Compose
using polyBindingModel

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("regression.jl")
include("synergy.jl")
include("translation.jl")
include("temporal.jl")

include("figures/figureCommon.jl")
include("figures/figureW.jl")
include("figures/figure1.jl")
include("figures/figure2.jl")
include("figures/figureB1.jl")
include("figures/figureDiseaseSpecific.jl")
include("figures/figureCelltypeSynergy.jl")

function figureAll()
    setGadflyTheme()

    figureB1()

    figure_Mmelanoma()
    figure_MITP()
    figure_Mblood()
    figure_Mbone()
    figure_MHIV()
    figure_MBcell()
    figure_Hblood()
    figure_Hbone()
    figure_Hspleen()
    figure_HITP()

    figure_Syn_ncMO()
    figure_Syn_cMO()
    figure_Syn_NKs()
    figure_Syn_Neu()
    figure_Syn_EO()
end

export polyfc, polyfc_ActV

end # module
