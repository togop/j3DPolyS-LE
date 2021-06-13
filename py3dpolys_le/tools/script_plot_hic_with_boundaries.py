# python channel libraries
import numpy as np
from matplotlib import pyplot as plt

# custom modules
import cooler_files as cool
import hic_arrays as arr


if __name__ == "__main__":

    # prepare invivo matrix
    filepath_invivo_hdf = "/media/cubix/D86E-6C50/hi-c/published/N2_chrX_hicpro_moushumi_20201002_5000.hdf5"
    filepath_invivo_cool = "/media/cubix/D86E-6C50/hi-c/published/N2_hicpro_moushumi_20201002_5000.cool"
    #array_invivo = arr.get_hic_array_from_hdf5(filepath_invivo_hdf, verbose=True)
    array_invivo = cool.get_cropped_single_chromosome_array_from_cool(filepath_invivo_cool)
    array_invivo = arr.normalize_invivo_array(array_invivo, verbose=True)
    #array_invivo = np.divide(array_invivo, 1.3)
    #arr.plot(array_invivo)

    # prepare mnl data
    filepath_mnl = "/media/cubix/D86E-6C50/simulations/boundary/boundary_mnl_bid/hic_003.hdf5"
    array_mnl = arr.get_hic_array_from_hdf5(filepath_mnl)
    print(array_mnl.shape, array_invivo.shape)
    shape_array = array_invivo.shape
    array_mnl = arr.resize(array_mnl, shape_array)
    array_mnl = arr.normalize_simulation_array(array_mnl, verbose=True)
    #array_mnl = np.multiply(array_mnl, 1.3)
    # arr.plot(array_mnl)

    # merge arrays and plot loops on it
    print("mnl array")
    a, b = arr.info(array_mnl, True)
    print("invivo array")
    a, b = arr.info(array_invivo, True)

    merged_array = arr.merge_two_hic_arrays(array_mnl, array_invivo)
    print("merged array")
    a, b = arr.info(merged_array, True)

    filepath_loops_boundaryfile_manual = "/media/cubix/D86E-6C50/hi-c/data_ce/boundary_sites/boundaries_chrX_manual_detection.csv"
    merged_array = arr.mark_loop_windows(merged_array, filepath_loops_boundaryfile_manual, verbose=False)
    arr.plot(merged_array, set_clim=(-3, 1))  # , "invivo_standard.png")
    print(merged_array, "\n\n\n")
    arr.save_with_pillow(merged_array)
    merged_array[np.where(merged_array == np.NINF)] = 0
    print(merged_array, "\n\n\n")
    merged_array = np.multiply(merged_array, -10)
    print(merged_array, "\n\n\n")


    # show generated plots
    plt.show()


