#!/usr/bin/env julia
# Verification script for j3DPolySLE package installation

println("Verifying j3DPolySLE package setup...")
println("=" ^ 50)

# Check Project.toml exists
if !isfile("Project.toml")
    error("Project.toml not found!")
end
println("✓ Project.toml found")

# Check main module file
if !isfile("src.jl/j3DPolySLE.jl")
    error("src/j3DPolySLE.jl not found!")
end
println("✓ Main module file found")

# Check all submodules exist
required_modules = [
    "Logging.jl",
    "JobRunner.jl",
    "HicAnalysis.jl",
    "PlotHic.jl",
    "PlotSimStats.jl",
    "HicConverters.jl",
    "Stats3dpolysLe.jl",
    "Runner3dpolysLe.jl"
]

for mod in required_modules
    if !isfile(joinpath("src.jl", mod))
        error("Missing module: $mod")
    end
    println("✓ $mod found")
end

# Check CLI scripts
cli_scripts = [
    "3dpolys_le_runner",
    "3dpolys_le_stats",
    "plot_hic",
    "plot_sim_stats",
    "hic_converters",
    "hdf5_to_cooler",
    "3dpolys_le"
]

for script in cli_scripts
    script_path = joinpath("bin", script)
    if !isfile(script_path)
        @warn "CLI script not found: $script_path"
    else
        println("✓ CLI script found: $script")
    end
end

# Try to load the package
println("\nTesting package loading...")
try
    using Pkg
    Pkg.activate(".")
    include("src.jl/j3DPolySLE.jl")
    println("✓ Package loads successfully")
catch e
    @error "Failed to load package" exception=(e, catch_backtrace())
    exit(1)
end

println("\n" * "=" ^ 50)
println("Package verification complete!")
println("\nTo install:")
println("  using Pkg")
println("  Pkg.develop(path=\".\")  # Development mode")
println("  # or")
println("  Pkg.add(path=\".\")      # Production mode")


