# python channel libraries
from matplotlib import pyplot as plt

# custom modules
import hic_arrays as arr


if __name__ == "__main__":  # this is to plot simulation only
    filepath_simulation = "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/hic_003.hdf5"

    array_simulation = arr.get_hic_array_from_hdf5(filepath_simulation)
    array_simulation = arr.normalize_simulation_array(array_simulation, verbose=True)

    print("write png")
    arr.plot_sim_for_figures(array_simulation,
                             "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/loading_mex_sym.png")
    print("done")


if __name__ == "fromcool__main__":  # this is to plot simulation only
    filepath_simulation = "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/rex/boundary_rex_bid_asym.2000.cool"

    array_simulation = arr.get_array_from_cool(filepath_simulation)
    array_simulation = arr.normalize_simulation_array(array_simulation, verbose=True)

    print("write png")
    arr.plot_sim_for_figures(array_simulation,
                             "/media/cubix/D86E-6C50/Thesis Shared/graphics/paper/rex/boundary_rex_bid_asym.png")
    print("done")