# import bioconda packages
import cooler
import h5py
import numpy as np
import os

# import custom modules
import hic_arrays as arr


def cool_to_arrayHDF5(filepath_input_cool, filename_output_hdf5):

    # read cooler data
    cooler_hic = cooler.Cooler(filepath_input_cool)
    shape_hic_matrix = cooler_hic.shape
    matrix_hic = cooler_hic.matrix(balance=False, sparse=False)

    # convert matrix to np.array
    array_hic = np.array(matrix_hic)

    # write numpy array to hdf5 file
    with h5py.File(filename_output_hdf5, 'w') as file_h5py:
        file_h5py.create_dataset("hic_map", shape_hic_matrix, dtype="i", data=array_hic)


def cooler_to_coarse_hdf5(filepath_input_cooler, filepath_generated_coarse_cooler, filepath_output_hdf, coarseninc_factor=2, verbose=False):

    if verbose:
        print("start coarsening cooler")
    cooler.coarsen_cooler(filepath_input_cooler, filepath_generated_coarse_cooler, coarseninc_factor, 500)
    if verbose:
        print("start convert coarse cooler to hic-array")
    cool_to_arrayHDF5(filepath_generated_coarse_cooler, filepath_output_hdf)


def get_single_chromosome_array_from_cool_alternative(filepath_input_cool="/mnt/imaging.data/gzala/published/N2_hicpro_moushumi_20201002_5000.cool",
                                          chromosome_to_crop="chrX"):

    cooler_hic = cooler.Cooler(filepath_input_cool)

    hic_array = cooler_hic.matrix(balance=True, sparse=True).fetch(chromosome_to_crop)

    return hic_array


def get_single_chromosome_array_from_cool(filepath_input_cool="/mnt/imaging.data/gzala/published/N2_hicpro_moushumi_20201002_5000.cool",
                                          chromosome_to_crop="chrX"):

    cooler_hic = cooler.Cooler(filepath_input_cool)

    hic_array = cooler_hic.matrix(balance=False, sparse=False).fetch(chromosome_to_crop)

    """ # hardbalance, first turn balance off in hic_array 
    weights = cooler_hic.bins().fetch('chrX')['weight']
    bias = weights.values

    hic_array.data = bias[hic_array.row] * bias[hic_array.col] * hic_array.data

    matrix_chrX_hard-balanced = hic_array.toarray()"""

    #  cooler_new_hic = hic_converters.hic_to_cool(hic_array, chromosome_to_crop, 5000, filepath_output_cool)

    """print(type(cooler_new_hic))

    matrix_uuh = hic_converters.cool_to_matrix(cooler_new_hic)

    print(type(matrix_uuh))

    mcool_new = hic_converters.matrix_to_mcool(matrix_uuh, 'chrX', 5000, [2, 4, 6, 8])

    print(type(mcool_new))

    file_mcool = open("N2_chrX_hicpro_moushumi_20201002_5000.mcool", "x")

    file_mcool.write(mcool_new)

    file_mcool.close()"""

    return hic_array


def get_single_chromosome_array_from_cool_alternative(filepath_input_cool="/mnt/imaging.data/gzala/published/N2_hicpro_moushumi_20201002_5000.cool",
                                          chromosome_to_crop="chrX"):

    cooler_hic = cooler.Cooler(filepath_input_cool)

    hic_array = cooler_hic.matrix(balance=True, sparse=True).fetch(chromosome_to_crop)

    return hic_array


def get_single_chromosome_array_from_cool_alternative(filepath_input_cool="/mnt/imaging.data/gzala/published/N2_hicpro_moushumi_20201002_5000.cool",
                                          chromosome_to_crop="chrX"):

    cooler_hic = cooler.Cooler(filepath_input_cool)

    hic_array = cooler_hic.matrix(balance=True, sparse=True).fetch(chromosome_to_crop)

    return hic_array


def get_single_chromosome_array_from_cool_alternative(filepath_input_cool="/mnt/imaging.data/gzala/published/N2_hicpro_moushumi_20201002_5000.cool",
                                          chromosome_to_crop="chrX"):

    cooler_hic = cooler.Cooler(filepath_input_cool)

    hic_array = cooler_hic.matrix(balance=True, sparse=False).fetch(chromosome_to_crop)

    return hic_array


def get_single_chromosome_balanced_array_from_cool(filepath_input_cool="/mnt/imaging.data/gzala/published/N2_hicpro_moushumi_20201002_5000.cool",
                                                   chromosome_to_crop="chrX"):

    cooler_hic = cooler.Cooler(filepath_input_cool)

    hic_matrix = cooler_hic.matrix(balance=False, sparse=True).fetch(chromosome_to_crop)

    weights = cooler_hic.bins().fetch(chromosome_to_crop)['weight']
    bias = weights.values

    hic_matrix.data = bias[hic_matrix.row] * bias[hic_matrix.col] * hic_matrix.data

    hic_array = hic_matrix.toarray()

    hic_array = arr.remove_nan_from_array(hic_array)

    return hic_array


def create_cropped_rebalanced_cool_from_whole_genome_cool(filepath_multiple_chromosome_cool, chromosome_to_crop, bin_size_of_original_cool):

    array_hic = get_single_chromosome_array_from_cool(filepath_multiple_chromosome_cool, chromosome_to_crop)

    arr.write_hic_array_to_cool(array_hic, chromosome_to_crop, bin_size_of_original_cool, filepath_multiple_chromosome_cool[:-4] + chromosome_to_crop + ".cool")

    os.system("cooler balance " + filepath_multiple_chromosome_cool[:-4] + chromosome_to_crop + ".cool")


if __name__ == "__main__":

    filepath_invivo_cool = "../../test/data/N2_hicpro_moushumi_20210618_2000.cool"

    create_cropped_rebalanced_cool_from_whole_genome_cool(filepath_invivo_cool, "chrX", 2000)


if __name__ == "__main__":
    #filepath_hdf5 = "/media/cubix/D86E-6C50/hi-c/scripts/N2_chrX_hicpro_moushumi_20201002_5000_coarse_test.hdf5"
    #filepath_cooler = "/media/cubix/D86E-6C50/hi-c/scripts/N2_chrX_hicpro_moushumi_20201002_5000_coarse_test.cool"

    filepath_hdf5 = "/media/cubix/D86E-6C50/hi-c/simulations/3DpolyS-LE/boundary/240_840k_sip_bd1insimstepalready/sim_sip_bd1insimstepalready.10000.hdf"
    filepath_cooler = "/media/cubix/D86E-6C50/hi-c/simulations/3DpolyS-LE/boundary/240_840k_sip_bd1insimstepalready/hic_003.hdf5.10000.cool"

    filepath_hardbalanced_invivo_coarse = "/media/cubix/D86E-6C50/hi-c/published/N2_chrX_hicpro_hard-balanced_moushumi_20201002_10000.cool"
    filepath_hardbalanced_invivo_coarse_hdf5 = "/media/cubix/D86E-6C50/hi-c/published/N2_chrX_hicpro_hard-balanced_moushumi_20201002_10000.hdf5"

    cool_to_arrayHDF5(filepath_hardbalanced_invivo_coarse, filepath_hardbalanced_invivo_coarse_hdf5)
