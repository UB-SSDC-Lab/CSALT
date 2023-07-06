
# Load dependencies
using FastGaussQuadrature
using StaticArrays
using ForwardDiff
using LinearAlgebra
using Match
using DelimitedFiles
using PyPlot

# Include source
include(joinpath(@__DIR__, "quadrature.jl"))
include(joinpath(@__DIR__, "eoms.jl"))
include(joinpath(@__DIR__, "analytic_partials.jl"))
include(joinpath(@__DIR__, "autodiff_partials.jl"))
include(joinpath(@__DIR__, "process_solution.jl"))
