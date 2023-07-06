
# This script tests that the analytic partial derivatives are correct

# Include source code
include(joinpath(@__DIR__,"..","include.jl"))

# Set relevant parameters
p       = 11359.07
f       = 0.7306
g       = 0.1
h       = 0.1
k       = 0.1
m       = 1000.0
t       = 0.0
δ       = 1.0
α1      = 1.0 / sqrt(3)
α2      = 1.0 / sqrt(3)
α3      = 1.0 / sqrt(3)
μ       = 3.986e5
tMax    = 0.5e-3
Isp     = 3000.0
g0      = 9.81e-3
τs,ws   = gausslegendre(10) 

# Compute partials
ddxdx       = averaged_reduced_state_dynamics_state_partial(
                    p,f,g,h,k,m,t,δ,α1,α2,α3,μ,tMax,Isp,g0,τs,ws)
ad_ddxdx    = ad_averaged_reduced_state_dynamics_state_partial(
                    p,f,g,h,k,m,t,δ,α1,α2,α3,μ,tMax,Isp,g0,τs,ws)
ddxdu       = averaged_reduced_state_dynamics_control_partial(
                    p,f,g,h,k,m,t,δ,α1,α2,α3,μ,tMax,Isp,g0,τs,ws)
ad_ddxdu    = ad_averaged_reduced_state_dynamics_control_partial(
                    p,f,g,h,k,m,t,δ,α1,α2,α3,μ,tMax,Isp,g0,τs,ws)
ddxdx_diff  = ddxdx .- ad_ddxdx
ddxdu_diff  = ddxdu .- ad_ddxdu

# Compute max diffs and print
ddxdx_maxdiff   = maximum(ddxdx_diff)
ddxdu_maxdiff   = maximum(ddxdu_diff)

println("Max state partial error:   " * string(ddxdx_maxdiff))
println("Max control partial error: " * string(ddxdu_maxdiff))