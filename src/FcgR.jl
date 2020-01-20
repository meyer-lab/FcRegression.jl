module FcgR
using LinearAlgebra
using ForwardDiff
using Calculus

include("fcBindingModel.jl")
include("dataHelpers.jl")
include("regression.jl")
include("synergy.jl")
include("translation.jl")

using Plots

include("figures/figure1.jl")
include("figures/figureB1.jl")
include("figures/figureB2.jl")

export calculateIsobologram, polyfc, polyfc_ActV

end # module
