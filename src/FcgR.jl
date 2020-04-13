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
include("systemsSerology.jl")

include("figures/figureW.jl")
include("figures/figure1.jl")
include("figures/figureB1.jl")
include("figures/figureB2.jl")

function figureAll()
    figureB1()
    figureB2()
    figureB3()
    figureB4()
    figureB5()
    figureB6()
    figureB11()
    figureB12()
    figureB13()
end

export polyfc, polyfc_ActV

end # module
