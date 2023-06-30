# The csaltApplicationWithJulia julia source code directory "jlsrc" is organized such that 
# all code corresponding to a specific CSALT problem be stored in its own directory, i.e.,
# "jlsrc/averaged_orbital_elements". This include file should only include each CSALT problems
# respective include file in its own directory.

# Set src base directory
const src_dir = @__DIR__

# Define functions for easy including of files
add_include_file(prob_src_dir)      = include(joinpath(src_dir, prob_src_dir, "include.jl"))
add_src_file(prob_src_dir, fname)   = include(joinpath(src_dir, prob_src_dir, fname))

# Include CSALT problem source code
add_include_file("averaged_orbital_elements")
