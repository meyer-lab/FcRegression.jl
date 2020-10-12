module FcRegression
using LinearAlgebra
import Distances
using Gadfly
using Compose
using polyBindingModel

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("regression.jl")
include("mixture.jl")
include("synergy.jl")
include("fitActI.jl")

include("figures/figureCommon.jl")
include("figures/figureW.jl")
include("figures/figure1.jl")
include("figures/figure2.jl")
include("figures/figure3.jl")
include("figures/figureB1.jl")
include("figures/figureS.jl")
include("figures/figureRR1.jl")
include("figures/figureDiseaseSpecific.jl")
include("figures/figureCelltypeSynergy.jl")
include("figures/figureDiseaseSynergy.jl")


function figureAll()
    setGadflyTheme()

    figureB1()
    figure2()
    figure3()
    plotMixtures()

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

    figure_MSyn_ncMO()
    figure_MSyn_cMO()
    figure_MSyn_NKs()
    figure_MSyn_Neu()
    figure_MSyn_EO()
    figure_HSyn_ncMO()
    figure_HSyn_cMO()
    figure_HSyn_NKs()
    figure_HSyn_Neu()
    figure_HSyn_EO()

    figure_MSyn_melanoma()
    figure_MSyn_ITP()

    figureRR1()
end

export polyfc, polyfc_ActV

end # module
