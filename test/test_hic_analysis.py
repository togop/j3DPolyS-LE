#! /usr/bin/env python

import os.path
import filecmp
import pytest
from py3dpolys_le import hic_analysis as ha


def test_chip_out_to_bedgraph():
    chip_out_file = './test/data/sim/Chip.out'
    chip_bedgraph_file = './out/Chip.bedGraph'
    ha.chip_out_to_bedgraph(chip_out_file, chip_bedgraph_file)
    expected = './test/expected/sim/Chip.bedGraph'
    assert os.path.exists(chip_bedgraph_file)
    assert filecmp.cmp(chip_bedgraph_file, expected, shallow=False)

