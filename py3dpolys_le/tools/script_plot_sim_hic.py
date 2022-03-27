# python channel libraries
from matplotlib import pyplot as plt

# custom modules
import hic_arrays as arr


if __name__ == "1__main__":  # this is to plot simulation only
    #filepath_simulation = "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/hic_003.hdf5"
    filepath_cool = "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/hic_003.hdf5"

    print("start")

    array_simulation = arr.get_hic_array_from_hdf5(filepath_cool)
    array_simulation = arr.normalize_simulation_array(array_simulation, verbose=True)

    print("write png")
    arr.plot_sim_for_figures(array_simulation,
                             "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/chrII_inserts/boundary_jimenez_chrII_tworexinsertions_sym.png")
    print("done")


if __name__ == "__main__":  # this is to plot simulation only
    filepath_cool = "/media/cubix/LaCie/ny_lab_dataset/chrII_3rex_insertions_17.cool"
    print("start")

    array_simulation = arr.get_array_from_cool(filepath_cool, "II")
    array_simulation = arr.normalize_simulation_array(array_simulation, verbose=True)

    print("write png")
    arr.plot_sim_for_figures(array_simulation,
                             "/media/cubix/LaCie/ny_lab_dataset/chrII_3rex_insertions_17.chrII.png")
    print("done")