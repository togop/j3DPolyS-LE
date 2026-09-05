#!/usr/bin/env julia
# Script to build HDF5.jl with correct system library paths
# Usage: julia --project=. build_hdf5.jl

using Pkg
Pkg.activate(".")

# Set HDF5_DIR to Homebrew installation
hdf5_prefix = readchomp(`brew --prefix hdf5`)
ENV["HDF5_DIR"] = hdf5_prefix

println("Setting HDF5_DIR to: $hdf5_prefix")
println("Building HDF5.jl...")

try
    Pkg.build("HDF5", verbose=true)
    println("✓ HDF5.jl built successfully!")
    
    # Test it
    using HDF5
    println("✓ HDF5.jl loaded successfully!")
catch e
    @error "Failed to build HDF5.jl" exception=(e, catch_backtrace())
    exit(1)
end

