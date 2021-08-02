import cooler
import hic_arrays as arr
import cooler_files as cool
import os


if __name__ == "__main__":
    filepath_mnl = "/media/cubix/D86E-6C50/simulations/boundary/asymetrical/boundary_manual_bid/hic_003.hdf5"
    filename = "boundary_manual_bid_asym.2000.cool"

    array_hic = arr.get_hic_array_from_hdf5(filepath_mnl)
    arr.write_hic_array_to_cool(array_hic, "chrX", 2000, filename)

    os.system("cooler zoomify -r 2000,4000,8000,16000,32000 " + filename)
