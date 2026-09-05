#!/bin/bash
# Benchmark test suite performance using time command

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║     TEST SUITE PERFORMANCE: JULIA vs PYTHON                    ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

cd "$(dirname "$0")"

# Julia Tests
echo "┌────────────────────────────────────────────────────────────────┐"
echo "│ JULIA TESTS (test/runtests.jl)                                │"
echo "└────────────────────────────────────────────────────────────────┘"
echo ""
echo "Running Julia tests 3 times..."
echo ""

for i in 1 2 3; do
    echo "Run $i:"
    /usr/bin/time -p julia --project=.. runtests.jl 2>&1 | grep -A3 "Test Summary" | head -4
    echo ""
done

echo ""
echo "┌────────────────────────────────────────────────────────────────┐"
echo "│ PYTHON TESTS (pytest)                                         │"
echo "└────────────────────────────────────────────────────────────────┘"
echo ""
echo "Note: Python tests need compiled binaries to run."
echo "      Showing collection time only..."
echo ""

source ../.venv/bin/activate

for i in 1 2 3; do
    echo "Run $i:"
    /usr/bin/time -p pytest test_3dpolys_le_runner.py test_hic_analysis.py --collect-only -q 2>&1 | tail -5
    echo ""
done

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                    DETAILED COMPARISON                         ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "To see detailed timings, run:"
echo ""
echo "  Julia:  time julia --project=. test/runtests.jl"
echo "  Python: time pytest test/ --collect-only"
echo ""
