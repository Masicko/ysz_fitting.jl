using Pkg
Pkg.add("Plots")
Pkg.add("VoronoiFVM")
Pkg.add("Revise")
Pkg.add("ForwardDiff")
# add
# 
# # try
# #         @eval using Revise
# #         # Turn on Revise s automatic evaluation behaviour
# #         Revise.async_steal_repl_backend()
# # catch err
# #         @warn "Could not load Revise."
# # end
#
# to  .julia/config/startup.jl 


Pkg.add("PyPlot")
Pkg.add("LsqFit") 
Pkg.add("PyPlot")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("LeastSquaresOptim")
Pkg.add("Optim")
Pkg.add("SparseArrays")
Pkg.add("ProgressMeter")
Pkg.add("LaTeXStrings")

# for NNLS
Pkg.add("MathProgBase")
# integration
Pkg.add("QuadGK")

# NEW
Pkg.add("BlackBoxOptim")
