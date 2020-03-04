module FcgR
using LinearAlgebra
using ForwardDiff
import Distances
using Gadfly
using GLM
using Compose
import Cairo

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("regression.jl")
include("synergy.jl")
include("translation.jl")


include("figures/figure1.jl")
include("figures/figureB1.jl")
include("figures/figureB2.jl")
include("figures/figureW.jl")

export polyfc, polyfc_ActV

end # module
