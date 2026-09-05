# Workaround script to load j3DPolySLE without precompilation
# This bypasses the precompilation error

using Pkg
Pkg.activate(".")

# Load the package directly without precompilation
include("src/j3DPolySLE.jl")
using .j3DPolySLE

println("✓ Package loaded successfully!")
println("You can now use:")
println("  j3DPolySLE.PlotHic.main()")
println("  j3DPolySLE.Runner3dpolysLe.main()")
println("  etc.")

