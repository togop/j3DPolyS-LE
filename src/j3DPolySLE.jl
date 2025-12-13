"""
    j3DPolySLE

    3D Polymer Simulations - Loop Extrusion model (Julia port)

    This package provides tools for simulating 3D polymer chromosome folding
    using loop extrusion factors (LEFs), boundary elements, and loading sites.

    # Usage

    ```julia
    using j3DPolySLE
    
    # Run simulations
    j3DPolySLE.Runner3dpolysLe.main()
    
    # Analyze statistics
    j3DPolySLE.Stats3dpolysLe.main()
    
    # Plot Hi-C data
    j3DPolySLE.PlotHic.main()
    ```

    See the README files for detailed documentation.
"""
module j3DPolySLE

# Export main modules
export Runner3dpolysLe, Stats3dpolysLe, HicAnalysis, JobRunner, PlotHic, PlotSimStats, HicConverters
export main

# Include all submodules
# Note: Logging.jl is a utility module but not exported to avoid conflicts
# All modules use the standard library Logging module directly
include("JobRunner.jl")
include("HicAnalysis.jl")
include("PlotHic.jl")
include("PlotSimStats.jl")
include("HicConverters.jl")
include("Stats3dpolysLe.jl")
include("Runner3dpolysLe.jl")

# Modules are automatically accessible as j3DPolySLE.ModuleName after include()
# The export statement makes them available when using j3DPolySLE

"""
    main()

Main entry point for the package. Calls the runner's main function.
"""
function main()
    Runner3dpolysLe.main()
end

end # module

