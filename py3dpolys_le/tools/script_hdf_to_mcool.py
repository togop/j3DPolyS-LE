import cooler
import hic_arrays as arr
import cooler_files as cool
import os


if __name__ == "__main__":
    filepath_mnl = "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/hic_003.hdf5"
    filename = "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/boundary_fusionXV_bid_sym.2000.cool"

    print("load hic array from .hdf file")
    array_hic = arr.get_hic_array_from_hdf5(filepath_mnl)
    print("write hic array to .cool file")
    arr.write_hic_array_to_cool(array_hic, "chrX", 2000, filename)
    print("done")
    #os.system("cooler zoomify -r 2000,4000,8000,16000,32000 " + filename)
