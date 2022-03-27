# python channel libraries
import numpy as np
from matplotlib import pyplot as plt

# custom modules
import cooler_files as cool
import hic_arrays as arr
import py3dpolys_le.hic_analysis as ha


if __name__ == "_1_main__":  # this is to plot simulation only
    filepath_invivo_cool = "../test/data/N2_hicpro_moushumi_20210618_2000.chrX.cool"
    filepath_mnl = "/tmp/hic_003.hdf5"
    """(chi2_lin, alpha_lin) = ha.compare_hic_chromosome(filepath_mnl, filepath_invivo_cool, chrs=['6', 'chrX', 'X'],
                                                      res=ha.RESOLUTION,
                                                      tads_boundary=None, norm=True,
                                                      chi2_mode=ha.CHI2_MODE_LINEAR, plot=False)
    array_invivo = cool.get_single_chromosome_array_from_cool(filepath_invivo_cool)
    array_invivo = arr.normalize_invivo_array(array_invivo, verbose=True)
    """
    array_mnl = arr.get_hic_array_from_hdf5(filepath_mnl)
    #shape_array = array_invivo.shape
    #array_mnl = arr.resize(array_mnl, shape_array)
    array_mnl = arr.normalize_simulation_array(array_mnl, verbose=True)

    #merged_array = arr.merge_two_hic_arrays(array_mnl, array_invivo)
    #merged_array = arr.mark_loop_windows(merged_array, filepath_loops_boundaryfile_manual, 2000, verbose=False)
    arr.plot_old(array_mnl, "boundary_mnl_bid_20210712.png")  # , "ecran_nlefx0.01.png")
    plt.show()


if __name__ == "2__main__":  # for plotting invivo only
    filepath_invivo_cool = "../test/data/N2_hicpro_moushumi_20210618_2000.chrX.cool"
    array_invivo = cool.get_single_chromosome_balanced_array_from_cool(filepath_invivo_cool)
    array_invivo = arr.normalize_invivo_array(array_invivo, alpha=16, verbose=True)
    arr.plot_old(array_invivo, "hybrid_b_rex_l_mex_20210709.png", (-16, 0))  # , "ecran_nlefx0.01.png")
    plt.show()


if __name__ == "1__main__":  # this is for testing invivo vs simulated

    # prepare invivo matrix
    #filepath_invivo_hdf = "/media/cubix/D86E-6C50/hi-c/published/N2_chrX_hicpro_moushumi_20201002_5000.hdf5"
    filepath_invivo_cool = "../test/data/N2_hicpro_moushumi_20210618_2000.chrX.cool"
    #filepath_mnl = "/media/cubix/D86E-6C50/simulations/boundary/boundary_mnl_bid/hic_003.hdf5"
    filepath_mnl = "/media/cubix/D86E-6C50/simulations/boundary/asymetrical/boundary_manual_bid/hic_003.hdf5"

    #array_invivo = arr.get_hic_array_from_hdf5(filepath_invivo_hdf, verbose=True)

    (chi2_lin, alpha_lin) = ha.compare_hic_chromosome(filepath_mnl, filepath_invivo_cool, chrs=['6', 'chrX', 'X'],
                                                      res=ha.RESOLUTION, tads_boundary=None, norm=True,
                                                      chi2_mode=ha.CHI2_MODE_LINEAR)
    print(alpha_lin)
    array_invivo = cool.get_single_chromosome_balanced_array_from_cool(filepath_invivo_cool)
    print("before nanremove", np.max(array_invivo), np.min(array_invivo))
    array_invivo = arr.remove_nan_from_array(array_invivo)
    print("after nanremove", np.max(array_invivo), np.min(array_invivo))
    array_invivo = array_invivo * alpha_lin
    print("after alphalin", np.max(array_invivo), np.min(array_invivo))
    array_invivo = np.log10(array_invivo)
    print("after log", np.max(array_invivo), np.min(array_invivo))

    #array_invivo = arr.normalize_invivo_array(array_invivo, alpha=alpha_lin, verbose=True)
    #array_invivo = np.divide(array_invivo, 1.3)
    #arr.plot(array_invivo)

    # prepare mnl data
    array_mnl = arr.get_hic_array_from_hdf5(filepath_mnl)
    print("\nmnl array before normalisation")
    a, b = arr.info_old(array_mnl, True)
    array_mnl = arr.remove_nan_from_array(array_mnl)
    array_mnl = np.log10(array_mnl)
    #print(array_mnl.shape, array_invivo.shape)
    #shape_array = array_invivo.shape
    #array_mnl = arr.resize(array_mnl, shape_array)
    #array_mnl = arr.normalize_simulation_array(array_mnl, verbose=True)
    #array_mnl = np.multiply(array_mnl, 1.3)
    # arr.plot(array_mnl)

    # merge arrays and plot loops on it
    print("\nmnl array")
    a, b = arr.info_old(array_mnl, True)
    print("invivo array")
    a, b = arr.info_old(array_invivo, True)

    merged_array = arr.merge_two_hic_arrays(array_mnl, array_invivo)

    print("merged array")
    a, b = arr.info_old(merged_array, True)

    filepath_simulation_loops_boundaryfile = "../py3dpolys_le/data/ce/boundary_sites/boundaries_chrX_manual_detection.csv"
    #merged_array = arr.mark_loop_windows(merged_array, filepath_loops_boundaryfile_manual, 2000, verbose=False)
    arr.plot_new(merged_array, "_boogiewoogie.png")

    # show generated plots
    plt.show()


if __name__ == "__main__":  # this is for testing invivo vs simulated

    # prepare invivo matrix
    filepath_invivo_cool = "../test/data/N2_hicpro_moushumi_20210618_2000.chrX.cool"
    filepath_mnl = "/media/cubix/D86E-6C50/simulations/boundary/asymetrical/boundary_manual_bid/hic_003.hdf5"
    filepath_simulation_loops_boundaryfile = "../py3dpolys_le/data/ce/boundary_sites/boundaries_chrX_manual_detection.csv"

    # set parameters
    alpha_lin = 16
    image_8bit_max = 255
    resolution = 2000

    # load arrays
    print("\n||| Load Hi-C Matrices |||")
    array_invivo = cool.get_single_chromosome_balanced_array_from_cool(filepath_invivo_cool)
    arr.info(array_invivo, "loaded experimental Hi-C from cool")

    array_simulation = arr.get_hic_array_from_hdf5(filepath_mnl)
    arr.info(array_simulation, "load simulation Hi-C from hdf")

    # transformations
    print("\n||| Transform Invivo Array |||")
    array_invivo = np.clip(array_invivo, 0, 0.005)
    arr.info(array_invivo, "afterclip")
    array_invivo = arr.transform_array_newmax(array_invivo, image_8bit_max)
    #median = arr.get_median_of_diagonal(array_invivo, diagonal_thickness=2)
    #array_invivo = arr.transform_array_newmax(array_invivo, image_8bit_max, median)
    arr.info(array_invivo, "after raise")
    lowest_value = arr.find_lowest_nonzero_value(array_invivo)
    array_invivo = np.clip(array_invivo, lowest_value, image_8bit_max)
    arr.info(array_invivo, "after second clip")
    array_invivo = arr.stretch_value_range(array_invivo)
    arr.info(array_invivo, "after stretch")
    #array_invivo = arr.contrast(array_invivo, 10)
    #arr.info(array_invivo, "after contrast")
    array_invivo = np.log10(array_invivo)
    arr.info(array_invivo, "after log")
    array_invivo = np.nan_to_num(array_invivo, posinf=image_8bit_max, neginf=0)
    arr.info(array_invivo, "after -inf removal")

    # plot
    #arr.plot_new(array_invivo, "_invivoonly_clip0.005_auto_contrast10.png", save=True)
    #arr.plot_new(array_invivo, "_invivoonly_clip0.01_auto.png", save=True)  #[1000:2000, 1000:2000], "invivo")
    #arr.histogram_new(array_invivo, verbose=True)

    print("\n||| Transform Simulation Array |||")
    array_simulation = arr.transform_array_newmax(array_simulation, image_8bit_max)
    arr.info(array_simulation, "after raise")
    array_simulation = np.log10(array_simulation)
    arr.info(array_simulation, "after log")
    array_simulation = np.nan_to_num(array_simulation, posinf=image_8bit_max, neginf=0)
    arr.info(array_simulation, "after ninf removal")

    # plot
    #arr.plot_new(array_simulation[1000:2000, 1000:2000], "simulation")

    # merge arrays and plot loops on it
    print("\n||| Merge Arrays |||")
    arr.info(array_simulation, "mnl array")
    arr.info(array_invivo, "invivo array")

    merged_array = arr.merge_two_hic_arrays(array_simulation, array_invivo)
    arr.info(merged_array, "merged array")

    print("\n||| Plot merged Array |||")
    merged_array = arr.mark_loop_windows(merged_array, filepath_simulation_loops_boundaryfile, resolution, verbose=False)
    arr.plot_new(merged_array, "_stretch_clip0.05_auto.marked.png", save=True)
    print("Done")

    # show generated plots
    # plt.show()
