#!/usr/bin/env python3
"""
Benchmark Python implementation performance
"""
import time
import numpy as np
from py3dpolys_le import hic_analysis as ha
from py3dpolys_le import _3dpolys_le_runner as runner

def benchmark_function(func, *args, iterations=100, **kwargs):
    """Benchmark a function with multiple iterations"""
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        times.append(end - start)

    return {
        'mean': np.mean(times),
        'std': np.std(times),
        'min': np.min(times),
        'max': np.max(times),
        'total': np.sum(times)
    }

def main():
    print("=" * 60)
    print("Python Performance Benchmarks")
    print("=" * 60)
    print()

    # Benchmark 1: remove_duplicates
    print("1. remove_duplicates (small list)")
    print("-" * 40)
    test_list = [1, 2, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8]
    stats = benchmark_function(ha.remove_duplicates, test_list, iterations=10000)
    print(f"   Mean time: {stats['mean']*1e6:.2f} µs")
    print(f"   Std dev:   {stats['std']*1e6:.2f} µs")
    print(f"   Min time:  {stats['min']*1e6:.2f} µs")
    print(f"   Max time:  {stats['max']*1e6:.2f} µs")
    print()

    # Benchmark 2: remove_duplicates with max
    print("2. remove_duplicates with max (small list)")
    print("-" * 40)
    test_list = list(range(100))
    stats = benchmark_function(ha.remove_duplicates, test_list, max=50, iterations=10000)
    print(f"   Mean time: {stats['mean']*1e6:.2f} µs")
    print(f"   Std dev:   {stats['std']*1e6:.2f} µs")
    print(f"   Min time:  {stats['min']*1e6:.2f} µs")
    print(f"   Max time:  {stats['max']*1e6:.2f} µs")
    print()

    # Benchmark 3: str2bool
    print("3. str2bool")
    print("-" * 40)
    stats = benchmark_function(runner.str2bool, "true", iterations=10000)
    print(f"   Mean time: {stats['mean']*1e6:.2f} µs")
    print(f"   Std dev:   {stats['std']*1e6:.2f} µs")
    print(f"   Min time:  {stats['min']*1e6:.2f} µs")
    print(f"   Max time:  {stats['max']*1e6:.2f} µs")
    print()

    # Benchmark 4: average_contact_prob (small matrix)
    print("4. average_contact_prob (10x10 matrix)")
    print("-" * 40)
    matrix_small = np.random.rand(10, 10)
    stats = benchmark_function(ha.average_contact_prob, matrix_small, 3, iterations=1000)
    print(f"   Mean time: {stats['mean']*1e6:.2f} µs")
    print(f"   Std dev:   {stats['std']*1e6:.2f} µs")
    print(f"   Min time:  {stats['min']*1e6:.2f} µs")
    print(f"   Max time:  {stats['max']*1e6:.2f} µs")
    print()

    # Benchmark 5: average_contact_prob (medium matrix)
    print("5. average_contact_prob (100x100 matrix)")
    print("-" * 40)
    matrix_medium = np.random.rand(100, 100)
    stats = benchmark_function(ha.average_contact_prob, matrix_medium, 30, iterations=100)
    print(f"   Mean time: {stats['mean']*1e3:.2f} ms")
    print(f"   Std dev:   {stats['std']*1e3:.2f} ms")
    print(f"   Min time:  {stats['min']*1e3:.2f} ms")
    print(f"   Max time:  {stats['max']*1e3:.2f} ms")
    print()

    # Benchmark 6: average_contact_prob (large matrix)
    print("6. average_contact_prob (1000x1000 matrix)")
    print("-" * 40)
    matrix_large = np.random.rand(1000, 1000)
    stats = benchmark_function(ha.average_contact_prob, matrix_large, 100, iterations=10)
    print(f"   Mean time: {stats['mean']*1e3:.2f} ms")
    print(f"   Std dev:   {stats['std']*1e3:.2f} ms")
    print(f"   Min time:  {stats['min']*1e3:.2f} ms")
    print(f"   Max time:  {stats['max']*1e3:.2f} ms")
    print()

    # Benchmark 7: get_decay_distribution (small)
    print("7. get_decay_distribution (10x10 matrix)")
    print("-" * 40)
    matrix_small = np.random.rand(10, 10)
    stats = benchmark_function(ha.get_decay_distribution, matrix_small,
                               confidence=0.0, min_dist=1, max_dist=8, iterations=100)
    print(f"   Mean time: {stats['mean']*1e3:.2f} ms")
    print(f"   Std dev:   {stats['std']*1e3:.2f} ms")
    print(f"   Min time:  {stats['min']*1e3:.2f} ms")
    print(f"   Max time:  {stats['max']*1e3:.2f} ms")
    print()

    # Benchmark 8: get_decay_distribution (medium)
    print("8. get_decay_distribution (100x100 matrix)")
    print("-" * 40)
    matrix_medium = np.random.rand(100, 100)
    stats = benchmark_function(ha.get_decay_distribution, matrix_medium,
                               confidence=0.0, min_dist=1, max_dist=80, iterations=10)
    print(f"   Mean time: {stats['mean']*1e3:.2f} ms")
    print(f"   Std dev:   {stats['std']*1e3:.2f} ms")
    print(f"   Min time:  {stats['min']*1e3:.2f} ms")
    print(f"   Max time:  {stats['max']*1e3:.2f} ms")
    print()

    # Benchmark 9: get_decay_distribution (large)
    print("9. get_decay_distribution (500x500 matrix)")
    print("-" * 40)
    matrix_large = np.random.rand(500, 500)
    stats = benchmark_function(ha.get_decay_distribution, matrix_large,
                               confidence=0.0, min_dist=1, max_dist=400, iterations=3)
    print(f"   Mean time: {stats['mean']:.2f} s")
    print(f"   Std dev:   {stats['std']:.2f} s")
    print(f"   Min time:  {stats['min']:.2f} s")
    print(f"   Max time:  {stats['max']:.2f} s")
    print()

    print("=" * 60)
    print("Python Benchmarks Complete")
    print("=" * 60)

if __name__ == "__main__":
    main()
