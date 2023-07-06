
# Load dependencies
using FastGaussQuadrature
using StaticArrays
using ForwardDiff

# Include source
include(joinpath(@__DIR__, "quadrature.jl"))
include(joinpath(@__DIR__, "eoms.jl"))
