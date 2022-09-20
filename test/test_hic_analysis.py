#! /usr/bin/env python

import os.path
import filecmp
import pytest
from py3dpolys_le import hic_analysis as ha


# @pytest.mark.trylast
# @pytest.mark.skip
def test_chip_out_to_bedgraph():
    chip_out_file = './test/data/sim/Chip.out'
    chip_bedgraph_file = './out/Chip.bedGraph'
    ha.chip_out_to_bedgraph(chip_out_file, bed_graph_file=chip_bedgraph_file)
    expected = './test/expected/sim/Chip.bedGraph'
    assert os.path.exists(chip_bedgraph_file)
    assert filecmp.cmp(chip_bedgraph_file, expected, shallow=False)


# @pytest.mark.skip
def test_chip_out_to_bedgraph_10kb():
    chip_out_file = './test/data/sim/Chip.out'
    chip_bedgraph_file = './out/Chip_10kb.bedGraph'
    ha.chip_out_to_bedgraph(chip_out_file, bed_graph_file=chip_bedgraph_file, chrom='chrZ', resolution=10000)
    expected = './test/expected/sim/Chip_10kb.bedGraph'
    assert os.path.exists(chip_bedgraph_file)
    assert filecmp.cmp(chip_bedgraph_file, expected, shallow=False)
