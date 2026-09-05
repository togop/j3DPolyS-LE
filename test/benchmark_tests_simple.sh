#!/bin/bash
# Simple benchmark of test suite performance

echo "════════════════════════════════════════════════════════════════"
echo "  TEST SUITE PERFORMANCE BENCHMARK: JULIA vs PYTHON"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Julia Tests
echo "Julia Tests (runtests.jl)"
echo "──────────────────────────────────────────────────────────────"
echo "Running Julia tests 3 times and measuring..."
echo ""

julia_times=""
for i in 1 2 3; do
    echo -n "  Run $i: "
    start=$(python3 -c "import time; print(time.perf_counter())")
    julia --project=.. runtests.jl > /dev/null 2>&1
    end=$(python3 -c "import time; print(time.perf_counter())")
    elapsed=$(python3 -c "print(round(($end - $start) * 1000, 1))")
    julia_times="$julia_times $elapsed"
    echo "${elapsed} ms"
done

julia_avg=$(python3 -c "times = [${julia_times}]; print(round(sum(times) / len(times), 1))")
echo ""
echo "  Average: ${julia_avg} ms"
echo ""

# Python Tests (collection only, since they need binaries)
echo "Python Tests (pytest collection)"
echo "──────────────────────────────────────────────────────────────"
echo "Running Python test collection 3 times..."
echo ""

source ../.venv/bin/activate

python_times=""
for i in 1 2 3; do
    echo -n "  Run $i: "
    start=$(python3 -c "import time; print(time.perf_counter())")
    pytest test_3dpolys_le_runner.py test_hic_analysis.py --collect-only -q > /dev/null 2>&1
    end=$(python3 -c "import time; print(time.perf_counter())")
    elapsed=$(python3 -c "print(round(($end - $start) * 1000, 1))")
    python_times="$python_times $elapsed"
    echo "${elapsed} ms"
done

python_avg=$(python3 -c "times = [${python_times}]; print(round(sum(times) / len(times), 1))")
echo ""
echo "  Average: ${python_avg} ms"
echo ""

# Summary
echo "════════════════════════════════════════════════════════════════"
echo "  SUMMARY"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "  Julia Tests:   ${julia_avg} ms (12 tests executed)"
echo "  Python Tests:  ${python_avg} ms (6 tests collected)"
echo ""

# Calculate which is faster
if (( $(echo "$julia_avg < $python_avg" | bc -l) )); then
    speedup=$(python3 -c "print(round(${python_avg} / ${julia_avg}, 2))")
    echo "  ✓ Julia is ${speedup}x faster"
else
    speedup=$(python3 -c "print(round(${julia_avg} / ${python_avg}, 2))")
    echo "  ✓ Python is ${speedup}x faster"
fi

echo ""
echo "Note: Python tests require compiled binaries to execute fully."
echo "      Comparison shows test framework overhead only."
echo ""
