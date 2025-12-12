"""
CoolerModule - Pure Julia implementation of cooler file format

This module provides a drop-in replacement for the Python cooler library,
specifically for the functions used in the 3DPolyS-LE project.
"""
module CoolerModule

include("Cooler.jl")
include("Zoomify.jl")

using .Cooler
using .Zoomify

# Re-export main functionality
export CoolerFile, CoolerMatrix
export create_cooler, zoomify_cooler
export chromnames, chromsizes, binsize, bins, matrix, fetch

# Python-compatible interface - removed to avoid infinite recursion
# Users should use CoolerFile directly

# Property-like access (for Python compatibility)
Base.getproperty(c::CoolerFile, s::Symbol) = begin
    if s === :chromnames
        return Cooler.chromnames(c)
    elseif s === :chromsizes
        return Cooler.chromsizes(c)
    elseif s === :binsize
        return Cooler.binsize(c)
    else
        getfield(c, s)
    end
end

# Allow calling bins() like a method
(c::CoolerFile)(; kwargs...) = Cooler.bins(c)

end # module
