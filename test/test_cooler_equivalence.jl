# Test suite for cooler functionality
using Test
using DataFrames
using HDF5

# Load the package
project_root = abspath(joinpath(@__DIR__, ".."))

# Patch doc! to handle abstract types gracefully
@eval Base.Docs begin
    doc!(::Type, ::Base.Docs.Binding, ::Base.Docs.DocStr) = nothing
    doc!(::Type, ::Any...) = nothing
end

include(joinpath(project_root, "src", "j3DPolySLE.jl"))
using .j3DPolySLE

@testset "Cooler Implementation Tests" begin

    @testset "Cooler Module Loading" begin
        @test isdefined(HicAnalysis, :CoolerModule)
        @test isdefined(HicAnalysis.CoolerModule, :CoolerFile)
        @test isdefined(HicAnalysis.CoolerModule, :create_cooler)
        @test isdefined(HicAnalysis.CoolerModule, :zoomify_cooler)
    end

    @testset "Create Simple Cooler File" begin
        # Create a simple test cooler file
        test_cool = joinpath(tempdir(), "test_cooler.cool")

        # Create bins
        bins = DataFrame()
        bins[!, :chrom] = ["chr1", "chr1", "chr1", "chr1"]
        bins[!, :start] = [0, 1000, 2000, 3000]
        bins[!, Symbol("end")] = [1000, 2000, 3000, 4000]

        # Create pixels (sparse format)
        pixels = Dict{String, Vector}()
        pixels["bin1_id"] = [0, 0, 1, 1, 2]
        pixels["bin2_id"] = [1, 2, 2, 3, 3]
        pixels["count"] = [10.0, 5.0, 15.0, 8.0, 12.0]

        # Create cooler file
        metadata = Dict{String, Any}(
            "bin-size" => 1000,
            "format" => "HDF5::Cooler",
            "format-version" => 3
        )

        HicAnalysis.CoolerModule.create_cooler(test_cool,
            bins=bins,
            pixels=pixels,
            metadata=metadata
        )

        @test isfile(test_cool)

        # Test reading the cooler file
        cooler = HicAnalysis.CoolerModule.CoolerFile("$test_cool::/")

        @test HicAnalysis.CoolerModule.binsize(cooler) == 1000
        @test HicAnalysis.CoolerModule.chromnames(cooler) == ["chr1"]

        sizes = HicAnalysis.CoolerModule.chromsizes(cooler)
        @test sizes["chr1"] == 4000

        # Test bins access
        bins_data = HicAnalysis.CoolerModule.bins(cooler)
        @test haskey(bins_data, "chrom")
        @test haskey(bins_data, "start")
        @test haskey(bins_data, "end")

        # Test matrix fetch
        mat_obj = HicAnalysis.CoolerModule.matrix(cooler, balance=false)
        mat = HicAnalysis.CoolerModule.Cooler.fetch(mat_obj, "chr1")

        @test size(mat) == (4, 4)
        @test mat[1, 2] == 10.0  # bin1_id=0, bin2_id=1
        @test mat[2, 1] == 10.0  # symmetric
        @test mat[1, 3] == 5.0   # bin1_id=0, bin2_id=2

        # Clean up
        rm(test_cool, force=true)
    end

    @testset "Utility Functions Match" begin
        # Test remove_duplicates
        @test HicAnalysis.remove_duplicates([1, 2, 2, 3, 3, 4]) == [1, 2, 3, 4]
        @test HicAnalysis.remove_duplicates([1, 2, 3, 4, 5], 3) == [1, 2, 3]

        # Test average_contact_prob
        test_matrix = Float64[
            1.0 2.0 3.0 4.0;
            2.0 1.0 2.0 3.0;
            3.0 2.0 1.0 2.0;
            4.0 3.0 2.0 1.0
        ]

        avg1, sd1 = HicAnalysis.average_contact_prob(test_matrix, 1)
        avg2, sd2 = HicAnalysis.average_contact_prob(test_matrix, 2)

        @test avg1 == 2.0
        @test avg2 == 3.0
    end

    @testset "Decay Distribution" begin
        test_hic = Float64[
            10.0 5.0 2.0 1.0 0.5;
            5.0 10.0 5.0 2.0 1.0;
            2.0 5.0 10.0 5.0 2.0;
            1.0 2.0 5.0 10.0 5.0;
            0.5 1.0 2.0 5.0 10.0
        ]

        dists, probs, _, _ = HicAnalysis.get_decay_distribution(test_hic, 0.0, 1, 4)

        @test dists == [1, 2, 3]
        @test probs == [5.0, 2.0, 1.0]
    end

    @testset "Chi2 Distance Range" begin
        # Test linear mode
        dist_range_linear = HicAnalysis.get_chi2_dist_range(HicAnalysis.CHI2_MODE_LINEAR, 10000)
        @test length(dist_range_linear) > 0
        @test all(d -> d > 0, dist_range_linear)
        @test issorted(dist_range_linear)

        # Test log mode
        dist_range_log = HicAnalysis.get_chi2_dist_range(HicAnalysis.CHI2_MODE_LOG, 10000)
        @test length(dist_range_log) > 0
        @test all(d -> d > 0, dist_range_log)
        @test issorted(dist_range_log)

        # At SIM_RESOLUTION, minimum should be CHI2_RANGE_START
        dist_range_sim = HicAnalysis.get_chi2_dist_range(HicAnalysis.CHI2_MODE_LINEAR, HicAnalysis.SIM_RESOLUTION)
        @test minimum(dist_range_sim) >= HicAnalysis.CHI2_RANGE_START
    end

    @testset "Fill Diagonal" begin
        test_mat = ones(5, 5)
        HicAnalysis.fill_diagonal!(test_mat, 0)

        for i in 1:5
            @test test_mat[i, i] == 0
        end

        # Non-diagonal elements should still be 1
        @test test_mat[1, 2] == 1.0
        @test test_mat[2, 1] == 1.0
    end
end

println("✅ All cooler implementation tests passed!")
