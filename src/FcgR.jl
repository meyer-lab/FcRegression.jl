module FcgR
using LinearAlgebra
using ForwardDiff
import Distances

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("regression.jl")
include("synergy.jl")
include("translation.jl")

include("figures/figureW.jl")
include("figures/figure1.jl")
include("figures/figureB1.jl")
include("figures/figureB2.jl")
include("figures/figureB3.jl")
include("figures/figureB4.jl")
include("figures/figureB5.jl")

export polyfc, polyfc_ActV

end # module
