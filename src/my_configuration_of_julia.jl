using Pkg
Pkg.add("Plots")
Pkg.add("VoronoiFVM")
Pkg.add("Revise")
Pkg.add("ForwardDiff")

################################################################
# add the following to >>  .julia/config/startup.jl
################################################################
#
# try
#         @eval using Revise
#         # Turn on Revise s automatic evaluation behaviour
#         Revise.async_steal_repl_backend()
# catch err
#         @warn "Could not load Revise."
# end
# 
################################################################

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
Pkg.add("ExtendableGrids")
Pkg.add("GridVisualize")
Pkg.add("DataFramesMeta")

# NEW
Pkg.add("BlackBoxOptim")
