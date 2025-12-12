# Test suite for j3DPolySLE
using Test
using Pkg

# Load the package directly since we're using src.jl/ instead of src/
# This bypasses Julia's package system which expects src/
project_root = abspath(joinpath(@__DIR__, ".."))
Pkg.activate(project_root)

# Patch doc! to handle abstract types gracefully - must be BEFORE loading modules
@eval Base.Docs begin
    doc!(::Type, ::Base.Docs.Binding, ::Base.Docs.DocStr) = nothing
    doc!(::Type, ::Any...) = nothing
end

# Load the package directly
include(joinpath(project_root, "src.jl", "j3DPolySLE.jl"))
using .j3DPolySLE

@testset "j3DPolySLE Tests" begin
    @testset "Module Loading" begin
        @test isdefined(j3DPolySLE, :Runner3dpolysLe)
        @test isdefined(j3DPolySLE, :Stats3dpolysLe)
        @test isdefined(j3DPolySLE, :HicAnalysis)
        @test isdefined(j3DPolySLE, :JobRunner)
        @test isdefined(j3DPolySLE, :PlotHic)
        @test isdefined(j3DPolySLE, :PlotSimStats)
        @test isdefined(j3DPolySLE, :HicConverters)
    end
    
    @testset "Constants" begin
        @test HicAnalysis.SIM_RESOLUTION == 2000
        @test HicAnalysis.RESOLUTION == 10000
        @test HicAnalysis.PLOT_FORMAT == "png"
    end
    
    @testset "Utility Functions" begin
        @test HicAnalysis.remove_duplicates([1, 2, 2, 3, 3, 4]) == [1, 2, 3, 4]
        @test HicAnalysis.remove_duplicates([1, 2, 3, 4, 5], 3) == [1, 2, 3]
    end
end


