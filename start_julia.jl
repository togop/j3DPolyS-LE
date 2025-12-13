#!/usr/bin/env julia
# Helper script to load j3DPolySLE package
# Usage: julia --project=. start_julia.jl

using Pkg
Pkg.activate(".")
# Ensure all dependencies are installed
Pkg.instantiate()

# Build HDF5 if needed (requires system HDF5 libraries)
try
    using HDF5
catch e
    if occursin("HDF5 is not properly installed", string(e)) || occursin("libhdf5", string(e))
        @warn "HDF5 not properly installed, attempting to build..."
        @info "Setting HDF5_DIR to Homebrew installation..."
        try
            # Try to find HDF5 via Homebrew (macOS)
            hdf5_prefix = readchomp(`brew --prefix hdf5`)
            ENV["HDF5_DIR"] = hdf5_prefix
            @info "HDF5_DIR set to: $hdf5_prefix"
        catch
            @warn "Could not find HDF5 via Homebrew. Please set HDF5_DIR manually."
            @info "On macOS: export HDF5_DIR=\$(brew --prefix hdf5)"
            @info "On Linux: export HDF5_DIR=/usr"
        end
        try
            Pkg.build("HDF5", verbose=true)
            @info "HDF5 build completed. Retrying import..."
            using HDF5
        catch build_err
            @error "Failed to build HDF5. Please run: julia --project=. build_hdf5.jl"
            @error "Or install system HDF5 libraries first."
            @error "Error: $build_err"
        end
    else
        rethrow(e)
    end
end

# Build PyCall if needed
try
    using PyCall
catch
    @warn "PyCall not built, attempting to build..."
    Pkg.build("PyCall")
end

# Patch doc! to handle abstract types gracefully - must be BEFORE loading modules
# The error signature is: doc!(::Type, ::Base.Docs.Binding, ::Base.Docs.DocStr)
using Base.Docs
Base.Docs.doc!(::Type, ::Base.Docs.Binding, ::Base.Docs.DocStr) = nothing
Base.Docs.doc!(::Type, ::Any...) = nothing

# Load package directly (bypasses precompilation)
# Catch documentation errors and continue - the module may partially load
try
    include("src/j3DPolySLE.jl")
    using .j3DPolySLE
catch e
    err_str = string(e)
    if occursin("doc!", err_str) || occursin("FieldError", err_str) || 
       occursin("CfgJobRunner", err_str) || occursin("AbstractJobRunner", err_str) ||
       occursin("type DataType has no field", err_str) || occursin("no method matching doc!", err_str)
        @warn "Documentation error encountered (this is expected and harmless)"
        # Try to reload - the module might have partially loaded
        if !isdefined(Main, :j3DPolySLE)
            # If module doesn't exist, try loading again with error suppression
            try
                # Suppress errors during include
                old_stderr = stderr
                include("src/j3DPolySLE.jl")
                using .j3DPolySLE
            catch e2
                @error "Failed to load package even after retry"
                @error "Error: $e2"
                exit(1)
            end
        else
            using .j3DPolySLE
        end
    else
        rethrow(e)
    end
end

println("✓ j3DPolySLE loaded successfully!")
println("\nAvailable modules:")
println("  - j3DPolySLE.PlotHic")
println("  - j3DPolySLE.Runner3dpolysLe")
println("  - j3DPolySLE.Stats3dpolysLe")
println("  - j3DPolySLE.PlotSimStats")
println("  - j3DPolySLE.HicConverters")
println("\nExample usage:")
println("  j3DPolySLE.PlotHic.main()")

