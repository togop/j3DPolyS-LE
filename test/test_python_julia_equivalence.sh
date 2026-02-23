#!/bin/bash
# Test script to verify Python and Julia implementations are equivalent

set -e

echo "========================================="
echo "Testing Python and Julia Equivalence"
echo "========================================="
echo ""

# Test 1: Module imports
echo "Test 1: Module Imports"
echo "----------------------"

# Python imports
source .venv/bin/activate
python3 -c "
from py3dpolys_le import hic_analysis
from py3dpolys_le import _3dpolys_le_runner
from py3dpolys_le import _3dpolys_le_stats
print('✓ Python modules imported successfully')
"

# Julia imports
julia --project=. -e "
using j3DPolySLE
println(\"✓ Julia modules imported successfully\")
"

echo ""

# Test 2: Utility Functions
echo "Test 2: Utility Functions (remove_duplicates, str2bool)"
echo "--------------------------------------------------------"

# Python tests
python3 << 'EOF'
from py3dpolys_le import hic_analysis as ha
from py3dpolys_le import _3dpolys_le_runner as runner

# remove_duplicates
assert ha.remove_duplicates([1, 2, 2, 3, 3, 4]) == [1, 2, 3, 4]
assert ha.remove_duplicates([1, 2, 3, 4, 5], max=3) == [1, 2, 3]

# str2bool
assert runner.str2bool("true") == True
assert runner.str2bool("false") == False

print("✓ Python utility functions work correctly")
EOF

# Julia tests
julia --project=. << 'EOF'
using j3DPolySLE

# remove_duplicates
@assert HicAnalysis.remove_duplicates([1, 2, 2, 3, 3, 4]) == [1, 2, 3, 4]
@assert HicAnalysis.remove_duplicates([1, 2, 3, 4, 5], 3) == [1, 2, 3]

# str2bool
@assert Runner3dpolysLe.str2bool("true") == true
@assert Runner3dpolysLe.str2bool("false") == false

println("✓ Julia utility functions work correctly")
EOF

echo ""

# Test 3: Average Contact Probability
echo "Test 3: Average Contact Probability Function"
echo "---------------------------------------------"

python3 << 'EOF'
from py3dpolys_le import hic_analysis as ha
import numpy as np

test_matrix = np.array([
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 1.0, 2.0, 3.0],
    [3.0, 2.0, 1.0, 2.0],
    [4.0, 3.0, 2.0, 1.0]
])

avg1, sd1 = ha.average_contact_prob(test_matrix, 1)
avg2, sd2 = ha.average_contact_prob(test_matrix, 2)

assert avg1 == 2.0 and sd1 == 0.0
assert avg2 == 3.0 and sd2 == 0.0

print(f"✓ Python: distance=1 gives avg={avg1}, distance=2 gives avg={avg2}")
EOF

julia --project=. << 'EOF'
using j3DPolySLE

test_matrix = Float64[
    1.0 2.0 3.0 4.0;
    2.0 1.0 2.0 3.0;
    3.0 2.0 1.0 2.0;
    4.0 3.0 2.0 1.0
]

avg1, sd1 = HicAnalysis.average_contact_prob(test_matrix, 1)
avg2, sd2 = HicAnalysis.average_contact_prob(test_matrix, 2)

@assert avg1 == 2.0 && sd1 == 0.0
@assert avg2 == 3.0 && sd2 == 0.0

println("✓ Julia: distance=1 gives avg=$avg1, distance=2 gives avg=$avg2")
EOF

echo ""

# Test 4: Decay Distribution
echo "Test 4: Decay Distribution Function"
echo "------------------------------------"

python3 << 'EOF'
from py3dpolys_le import hic_analysis as ha
import numpy as np

test_hic = np.array([
    [10.0, 5.0, 2.0, 1.0, 0.5],
    [5.0, 10.0, 5.0, 2.0, 1.0],
    [2.0, 5.0, 10.0, 5.0, 2.0],
    [1.0, 2.0, 5.0, 10.0, 5.0],
    [0.5, 1.0, 2.0, 5.0, 10.0]
])

dists, probs, _, _ = ha.get_decay_distribution(test_hic, confidence=0.0, min_dist=1, max_dist=4)

assert list(dists) == [1, 2, 3]
assert [float(p) for p in probs] == [5.0, 2.0, 1.0]

print(f"✓ Python: distances={list(dists)}, probs={[float(p) for p in probs]}")
EOF

julia --project=. << 'EOF'
using j3DPolySLE

test_hic = Float64[
    10.0 5.0 2.0 1.0 0.5;
    5.0 10.0 5.0 2.0 1.0;
    2.0 5.0 10.0 5.0 2.0;
    1.0 2.0 5.0 10.0 5.0;
    0.5 1.0 2.0 5.0 10.0
]

dists, probs, _, _ = HicAnalysis.get_decay_distribution(test_hic, 0.0, 1, 4)

@assert dists == [1, 2, 3]
@assert probs == [5.0, 2.0, 1.0]

println("✓ Julia: distances=$dists, probs=$probs")
EOF

echo ""
echo "========================================="
echo "✅ ALL TESTS PASSED!"
echo "Python and Julia implementations are equivalent"
echo "========================================="
