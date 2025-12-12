#!/usr/bin/env julia
"""
Benchmark Julia implementation performance
"""

using Statistics
using Random

# Load the package
project_root = abspath(joinpath(@__DIR__, ".."))
include(joinpath(project_root, "src", "j3DPolySLE.jl"))
using .j3DPolySLE

function benchmark_function(func, args...; iterations=100, kwargs...)
    """Benchmark a function with multiple iterations"""
    times = Float64[]

    # Warm-up run
    func(args...; kwargs...)

    for _ in 1:iterations
        start = time_ns()
        result = func(args...; kwargs...)
        elapsed = time_ns() - start
        push!(times, elapsed / 1e9)  # Convert to seconds
    end

    return Dict(
        "mean" => mean(times),
        "std" => std(times),
        "min" => minimum(times),
        "max" => maximum(times),
        "total" => sum(times)
    )
end

function main()
    println("=" ^ 60)
    println("Julia Performance Benchmarks")
    println("=" ^ 60)
    println()

    # Benchmark 1: remove_duplicates
    println("1. remove_duplicates (small list)")
    println("-" ^ 40)
    test_list = [1, 2, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8]
    stats = benchmark_function(HicAnalysis.remove_duplicates, test_list; iterations=10000)
    println("   Mean time: $(round(stats["mean"]*1e6, digits=2)) µs")
    println("   Std dev:   $(round(stats["std"]*1e6, digits=2)) µs")
    println("   Min time:  $(round(stats["min"]*1e6, digits=2)) µs")
    println("   Max time:  $(round(stats["max"]*1e6, digits=2)) µs")
    println()

    # Benchmark 2: remove_duplicates with max
    println("2. remove_duplicates with max (small list)")
    println("-" ^ 40)
    test_list = collect(0:99)
    stats = benchmark_function(HicAnalysis.remove_duplicates, test_list, 50; iterations=10000)
    println("   Mean time: $(round(stats["mean"]*1e6, digits=2)) µs")
    println("   Std dev:   $(round(stats["std"]*1e6, digits=2)) µs")
    println("   Min time:  $(round(stats["min"]*1e6, digits=2)) µs")
    println("   Max time:  $(round(stats["max"]*1e6, digits=2)) µs")
    println()

    # Benchmark 3: str2bool
    println("3. str2bool")
    println("-" ^ 40)
    stats = benchmark_function(Runner3dpolysLe.str2bool, "true"; iterations=10000)
    println("   Mean time: $(round(stats["mean"]*1e6, digits=2)) µs")
    println("   Std dev:   $(round(stats["std"]*1e6, digits=2)) µs")
    println("   Min time:  $(round(stats["min"]*1e6, digits=2)) µs")
    println("   Max time:  $(round(stats["max"]*1e6, digits=2)) µs")
    println()

    # Benchmark 4: average_contact_prob (small matrix)
    println("4. average_contact_prob (10x10 matrix)")
    println("-" ^ 40)
    Random.seed!(42)
    matrix_small = rand(10, 10)
    stats = benchmark_function(HicAnalysis.average_contact_prob, matrix_small, 3; iterations=1000)
    println("   Mean time: $(round(stats["mean"]*1e6, digits=2)) µs")
    println("   Std dev:   $(round(stats["std"]*1e6, digits=2)) µs")
    println("   Min time:  $(round(stats["min"]*1e6, digits=2)) µs")
    println("   Max time:  $(round(stats["max"]*1e6, digits=2)) µs")
    println()

    # Benchmark 5: average_contact_prob (medium matrix)
    println("5. average_contact_prob (100x100 matrix)")
    println("-" ^ 40)
    Random.seed!(42)
    matrix_medium = rand(100, 100)
    stats = benchmark_function(HicAnalysis.average_contact_prob, matrix_medium, 30; iterations=100)
    println("   Mean time: $(round(stats["mean"]*1e3, digits=3)) ms")
    println("   Std dev:   $(round(stats["std"]*1e3, digits=3)) ms")
    println("   Min time:  $(round(stats["min"]*1e3, digits=3)) ms")
    println("   Max time:  $(round(stats["max"]*1e3, digits=3)) ms")
    println()

    # Benchmark 6: average_contact_prob (large matrix)
    println("6. average_contact_prob (1000x1000 matrix)")
    println("-" ^ 40)
    Random.seed!(42)
    matrix_large = rand(1000, 1000)
    stats = benchmark_function(HicAnalysis.average_contact_prob, matrix_large, 100; iterations=10)
    println("   Mean time: $(round(stats["mean"]*1e3, digits=3)) ms")
    println("   Std dev:   $(round(stats["std"]*1e3, digits=3)) ms")
    println("   Min time:  $(round(stats["min"]*1e3, digits=3)) ms")
    println("   Max time:  $(round(stats["max"]*1e3, digits=3)) ms")
    println()

    # Benchmark 7: get_decay_distribution (small)
    println("7. get_decay_distribution (10x10 matrix)")
    println("-" ^ 40)
    Random.seed!(42)
    matrix_small = rand(10, 10)
    stats = benchmark_function(HicAnalysis.get_decay_distribution, matrix_small,
                               0.0, 1, 8; iterations=100)
    println("   Mean time: $(round(stats["mean"]*1e3, digits=2)) ms")
    println("   Std dev:   $(round(stats["std"]*1e3, digits=2)) ms")
    println("   Min time:  $(round(stats["min"]*1e3, digits=2)) ms")
    println("   Max time:  $(round(stats["max"]*1e3, digits=2)) ms")
    println()

    # Benchmark 8: get_decay_distribution (medium)
    println("8. get_decay_distribution (100x100 matrix)")
    println("-" ^ 40)
    Random.seed!(42)
    matrix_medium = rand(100, 100)
    stats = benchmark_function(HicAnalysis.get_decay_distribution, matrix_medium,
                               0.0, 1, 80; iterations=10)
    println("   Mean time: $(round(stats["mean"]*1e3, digits=2)) ms")
    println("   Std dev:   $(round(stats["std"]*1e3, digits=2)) ms")
    println("   Min time:  $(round(stats["min"]*1e3, digits=2)) ms")
    println("   Max time:  $(round(stats["max"]*1e3, digits=2)) ms")
    println()

    # Benchmark 9: get_decay_distribution (large)
    println("9. get_decay_distribution (500x500 matrix)")
    println("-" ^ 40)
    Random.seed!(42)
    matrix_large = rand(500, 500)
    stats = benchmark_function(HicAnalysis.get_decay_distribution, matrix_large,
                               0.0, 1, 400; iterations=3)
    println("   Mean time: $(round(stats["mean"], digits=4)) s")
    println("   Std dev:   $(round(stats["std"], digits=4)) s")
    println("   Min time:  $(round(stats["min"], digits=4)) s")
    println("   Max time:  $(round(stats["max"], digits=4)) s")
    println()

    println("=" ^ 60)
    println("Julia Benchmarks Complete")
    println("=" ^ 60)
end

main()
