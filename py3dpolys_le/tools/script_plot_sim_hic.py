# python channel libraries
import numpy as np
from matplotlib import pyplot as plt

# custom modules
import cooler_files as cool
import hic_arrays as arr
import py3dpolys_le.hic_analysis as ha


if __name__ == "__main__":  # this is to plot simulation only
    filepath_invivo_cool = "../test/data/N2_hicpro_moushumi_20210618_2000.chrX.cool"
    filepath_mnl = "/tmp/hic_003.hdf5"

    array_mnl = arr.get_hic_array_from_hdf5(filepath_mnl)
    #shape_array = array_invivo.shape
    #array_mnl = arr.resize(array_mnl, shape_array)
    array_mnl = arr.normalize_simulation_array(array_mnl, verbose=True)

    #merged_array = arr.merge_two_hic_arrays(array_mnl, array_invivo)
    #merged_array = arr.mark_loop_windows(merged_array, filepath_loops_boundaryfile_manual, 2000, verbose=False)
    arr.plot_old(array_mnl, "hybrid_l_mex_b_rex_bid_sym.png")  # , "ecran_nlefx0.01.png")
    plt.show()
