# The csaltApplicationWithJulia julia source code directory "jlsrc" is organized such that 
# all code corresponding to a specific CSALT problem be stored in its own directory, i.e.,
# "jlsrc/AveragedOrbitalElements". 

# First, we need to activate and instantiate the Julia CSALT environment
using Pkg;
Pkg.activate(@__DIR__);
Pkg.instantiate();

# Next, include project files
include(joinpath(@__DIR__, "AveragedOrbitalElements", "include.jl"))