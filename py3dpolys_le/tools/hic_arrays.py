import h5py
import cooler
import numpy as np
import pandas as pd
from PIL import Image
from matplotlib import pyplot as plt

import loop_files as loops


def get_hic_array_from_hdf5(filepath_array_hdf5, verbose=False):

    if verbose:
        print("load hic-array from hdf: ", filepath_array_hdf5)

    # show hic map
    with h5py.File(filepath_array_hdf5, 'r') as file_h5py:

        # do some accessor magic to read data
        key_hic_map__h5py = list(file_h5py.keys())[0]
        list_hic = list(file_h5py[key_hic_map__h5py])
        hic_array = np.array(list_hic)

    return hic_array


def info(array_hic, verbose=False):
    # get max/min value in method array
    val_max = np.max(array_hic)
    val_min = np.min(array_hic)
    shape = array_hic.shape
    if verbose:
        print("Hi_C:: max", val_max, " / min", val_min, "| Shape: ", shape)

    return val_max, val_min


def plot(hic_array, filename_output=None, set_clim=(-2.75, 0)):

    val_max, val_min = info(hic_array)

    # show hic image
    fig = plt.figure()

    plt.imshow(hic_array, interpolation='nearest', cmap="hot_r")
    plt.clim(set_clim)

    plt.colorbar()

    if filename_output is not None:
        fig.savefig(filename_output, dpi=2000)


def get_array_from_cool(filepath_input_cool):

    cooler_hic = cooler.Cooler(filepath_input_cool)

    matrix = cooler_hic.matrix(balance=False, sparse=False).fetch('chrX')
    print(np.shape(matrix))

    print(matrix.data)

    #hic_array = matrix.toarray()

    return matrix


def get_balanced_array_from_cool(filepath_input_cool):

    cooler_hic = cooler.Cooler(filepath_input_cool)

    # get sparse matrix
    matrix = cooler_hic.matrix(balance=False, sparse=True).fetch('chrX')

    # get bias from weights
    weights = cooler_hic.bins().fetch('chrX')['weight']
    bias = weights.values

    # do the balancing
    matrix.data = bias[matrix.row] * bias[matrix.col] * matrix.data

    # convert sparse matrix to array
    matrix_chrX_hardbalanced = matrix.toarray()

    return matrix_chrX_hardbalanced


def write_hic_array_to_hdf(hic_array, filepath_output_hdf, verbose=False):

    shape_hic_matrix = np.shape(hic_array)

    if verbose:
        print("start writing process")
        print("write .hdf with shape", shape_hic_matrix, " to ", filepath_output_hdf)

    with h5py.File(filepath_output_hdf, 'w') as file_h5py:
        file_h5py.create_dataset("hic_map", shape_hic_matrix, dtype="f", data=hic_array)


def remove_zeros_from_hic_array(hic_array):

    indices_zeroes = np.where(np.isnan(hic_array))

    print("indices\n", indices_zeroes)

    hic_array[indices_zeroes] = 0

    return hic_array


def write_hic_array_to_cool(hic, chr, resolution, cool_file, verbose=False):

    if verbose:
        print("make cool from hic-array")

    # build the cooler fields
    N = hic.shape[0]

    bins_index = [[chr, i * resolution, i * resolution + resolution] for i in range(N)]
    bins = pd.DataFrame(data=bins_index, columns=['chrom', 'start', 'end'])  # , dtype=np.dtype([('','','')]))

    pixels_bin1_id = []
    pixels_bin2_id = []
    pixels_count = []

    tot_iter = (N - 1) * N / 2
    iter = 0
    for bin1_id in range(N - 1):
        for bin2_id in range(bin1_id + 1, N):
            iter += 1
            progress = (iter / tot_iter) * 100


            count = hic[bin1_id, bin2_id]

            if count != 0:
                # pixels_pd = pixels_pd.append({'bin1_id': np.int64(bin1_id), 'bin2_id': np.int64(bin2_id), 'count': count}, ignore_index=True)
                pixels_bin1_id.append(np.int64(bin1_id))
                pixels_bin2_id.append(np.int64(bin2_id))
                pixels_count.append(count)

    pixels_dic = {'bin1_id': pixels_bin1_id, 'bin2_id': pixels_bin2_id, 'count': pixels_count}
    metadata = {'format': 'HDF5::Cooler',
                'format-version': '0.8.6',
                'bin-type': 'fixed',
                'bin-size': resolution,
                'storage-mode': 'symmetric-upper',
                'genome-assembly': 'ce11',
                'generated-by': "gabriels scripts" + '-' + "0.01",
                # 'creation-date': datetime.date.today()
                }

    count_dtypes = {'count': 'float64'}
    cooler.create_cooler(cool_file, bins=bins, pixels=pixels_dic, dtypes=count_dtypes, ordered=True, metadata=metadata)
    return cool_file
    # problem with showing .cool file in higlass but with .mcool it works


def normalize_invivo_array(hic_array, verbose=False):

    import py3dpolys_le.hic_analysis as ha
    # get alpha_lin
    sim_hic_file = "/media/cubix/D86E-6C50/simulations/boundary/boundary_mnl_bid/hic_003.hdf5"
    cmp_hic_file = "/media/cubix/D86E-6C50/hi-c/published/N2_hicpro_moushumi_20201002_5000.hdf"
    (chi2_lin, alpha_lin) = ha.compare_hic_chromosome(sim_hic_file, cmp_hic_file, chrs=['6', 'chrX', 'X'],
                                                      res=ha.RESOLUTION,
                                                      tads_boundary=None, norm=True,
                                                      chi2_mode=ha.CHI2_MODE_LINEAR)  # , plot=True)
    print(alpha_lin)

    if verbose:
        print("\nstart normalizing invivo hic array with proprieties:")
        unused_var_1, unused_var_2 = info(hic_array, True)
        print("remove diagonal from hic array")

    # fill the empty middle lane with max value in array
    np.fill_diagonal(hic_array, 0)

    # get data from array
    if verbose:
        val_max, val_min = info(hic_array, True)
    else:
        val_max, val_min = info(hic_array)

    # divide matrix so it has val_max = 1
    if verbose:
        print("divide array by factor", val_max)
    #hic_array = np.divide(hic_array, val_max)

    hic_array = np.divide(hic_array, alpha_lin * 10)

    # log transform matrix
    if verbose:
        unused_var_1, unused_var_2 = info(hic_array, True)
        print("apply log10 transform on array")
    hic_array = np.log10(hic_array)

    if verbose:
        print("\ndone normalizing invivo hic array")
        unused_var_1, unused_var_2 = info(hic_array, True)
        print("-------------------")

    return hic_array


def normalize_simulation_array(hic_array, verbose=False):

    if verbose:
        print("\nstart normalizing simulation hic array with proprieties:")
        unused_var_1, unused_var_2 = info(hic_array, True)
        print("remove diagonal from hic array")

    # fill the empty middle lane with max value in array
    np.fill_diagonal(hic_array, 0)

    # get data from array
    if verbose:
        val_max, val_min = info(hic_array, True)
    else:
        val_max, val_min = info(hic_array)

    # divide matrix so it has val_max = 1
    if verbose:
        print("divide array by factor", val_max)
    hic_array = np.divide(hic_array, val_max)

    # log transform matrix
    if verbose:
        unused_var_1, unused_var_2 = info(hic_array, True)
        print("apply log10 transform on array")
    hic_array = np.log10(hic_array)

    if verbose:
        print("\ndone normalizing simulation hic array")
        unused_var_1, unused_var_2 = info(hic_array, True)
        print("-------------------")

    return hic_array


def merge_two_hic_arrays(hic_array_upper, hic_array_lower):

    # delete bottom half of hic-array
    for index_row in range(len(hic_array_upper)):
        hic_array_upper[index_row][:index_row] = hic_array_lower[index_row][:index_row]

    # get data from array
    val_max, val_min = info(hic_array_upper, True)
    val_max, val_min = info(hic_array_lower, True)


    #print(np.where(hic_array == val_max))

    # fill the empty middle lane with max value in array
    np.fill_diagonal(hic_array_upper, val_max)

    return hic_array_upper


def resize(hic_array, shape_array):

    hic_array = np.multiply(hic_array, 1000)
    image_hic = Image.fromarray(hic_array)
    image_hic = image_hic.convert("L")
    image_resized = image_hic.resize(shape_array)
    hic_array = np.array(image_resized)

    return hic_array


def save_with_pillow(hic_array):

    image_hic = Image.fromarray(hic_array)
    image_hic.save("hic_array.tiff")


def mark_loop_windows(array_hic, filepath_boundaryfile_csv, verbose=False):

    if verbose:
        print("\nstart marking loop windows on hic array")

    # get data from array
    if verbose:
        print("input array data:")
        val_max, val_min = info(array_hic, True)
        shape_array_hic = np.shape(array_hic)
        index_max_array_hic = shape_array_hic[0]

    else:
        val_max, val_min = info(array_hic, False)
        shape_array_hic = np.shape(array_hic)
        index_max_array_hic = shape_array_hic[0]

    # mark anchors and loops
    thickness_line = 0
    wideness_window = 5

    # get loop coordinates
    dict_loops = loops.get_dict_from_boundaryfile(filepath_boundaryfile_csv, resolution=5000)
    if verbose:
        print("\n-------------------------\nLoops found in input:")

    ###########################################
    # draw loop-windows on hic array
    list_window_margins = []  # this list is used to make line with cutouts for windows
    for key in dict_loops.keys():
        list_bin_indices = sorted(dict_loops[key], reverse=True)

        if verbose:
            print("loop#", key, ": y =", list_bin_indices[0], "/ x = ", list_bin_indices[1])

        loop_y_pos = list_bin_indices[0]
        loop_x_pos = list_bin_indices[1]

        # window coordinates
        window_min_y = loop_y_pos - wideness_window
        window_max_y = loop_y_pos + wideness_window
        # print(window_max_y, window_min_y)

        window_min_x = loop_x_pos - wideness_window
        window_max_x = loop_x_pos + wideness_window
        # print(window_max_x, window_min_x)

        # set to max_index if window would be out of index
        if window_min_y >= index_max_array_hic:
            window_min_y = index_max_array_hic - 2

        if window_max_y >= index_max_array_hic:
            window_max_y = index_max_array_hic - 2

        if window_min_x >= index_max_array_hic:
            window_min_y = index_max_array_hic - 2

        if window_max_x >= index_max_array_hic:
            window_max_x = index_max_array_hic - 2

        # write window margins to list, y on index 0, x on index 1
        window_margins = np.array([[window_min_y, window_max_y], [window_min_x, window_max_x]])
        list_window_margins.append(window_margins)

        # plot top window line
        array_hic[window_min_y, window_min_x:window_max_x] = val_max

        # plot bot window line
        array_hic[window_max_y, window_min_x:window_max_x + 1] = val_max  # the 1 here is needed for nice squares

        # plot left window line
        array_hic[window_min_y:window_max_y, window_min_x] = val_max

        # plot left window line
        array_hic[window_min_y:window_max_y, window_max_x] = val_max

    if verbose:
        print("\n-------------------------\nlist window margins")
        counter = 0
        for window_margins in list_window_margins:
            print("window", counter, ":: height =", window_margins[0], "; width =", window_margins[1])
            counter += 1
        counter = 0  # safetyreset

    ###########################################
    # get subtraction intervals as dict of sets for vertical lines
    dict_vertical_line_cutouts = {}
    dict_horizontal_line_cutouts = {}
    for key in dict_loops.keys():
        list_bin_indices = sorted(dict_loops[key], reverse=True)

        if verbose:
            pass  # print("loop#", key, ": y =", list_bin_indices[0], "/ x = ", list_bin_indices[1])

        loop_y_pos = list_bin_indices[0]
        loop_x_pos = list_bin_indices[1]
        for window_margins in list_window_margins:
            margin_y = range(window_margins[0][0], window_margins[0][1])
            margin_x = range(window_margins[1][0], window_margins[1][1])
            if loop_x_pos in margin_x:

                # make sets in the dictionary for vertical lines
                if loop_x_pos in dict_vertical_line_cutouts.keys():
                    dict_vertical_line_cutouts[loop_x_pos].add(margin_y)

                else:
                    dict_vertical_line_cutouts[loop_x_pos] = {margin_y}

            if loop_y_pos in margin_y:
                # make sets in the dictionary
                if loop_y_pos in dict_horizontal_line_cutouts.keys():
                    dict_horizontal_line_cutouts[loop_y_pos].add(margin_x)

                else:
                    dict_horizontal_line_cutouts[loop_y_pos] = {margin_x}

    if verbose:
        print("-------------------------\nDict Substraction Intervals for each X-axis position:")
        for key in dict_vertical_line_cutouts.keys():
            print("vertical line on bin", key, "needs following cutouts:\n", dict_vertical_line_cutouts[key])

    ###########################################
    # cut vertical lines that cross a window and mark each line fragment on the array
    if verbose:
        print("\n-------------------------\nVertical lines")

    for key in dict_loops.keys():
        list_bin_indices = sorted(dict_loops[key], reverse=True)

        if verbose:
            print("\nVertical Line for loop #", key, "has starts and stops at y-axis bin indices:", list_bin_indices)

        v_line_min = list_bin_indices[1]
        v_line_max = list_bin_indices[0]

        v_line_y_position = v_line_min

        if verbose:
            print("y min: ", v_line_min, "y max: ", v_line_max)

        set_subtraction_intervals = dict_vertical_line_cutouts[v_line_y_position]
        if verbose:
            print("list subrtraction intervals vor vertical line on X position =", v_line_y_position)
            print(set_subtraction_intervals)

        # start getting the lines to plot
        list_unsorted_subtraction_intervals = list(set_subtraction_intervals)
        list_subtraction_intervals = sorted(list_unsorted_subtraction_intervals, key=lambda r: r.start)
        n_subtraction_intervals = len(list_subtraction_intervals)

        for index_interval in range(len(list_subtraction_intervals)):

            interval = list_subtraction_intervals[index_interval]

            if verbose:
                print("removing inverval", interval[0], "-", interval[-1], "from vertical line at Xpos =", v_line_y_position)

            # plot first vertical line fragment of line
            if index_interval == 0:

                start_horizontal_line = v_line_min
                end_horizontal_line_fragment = interval[0]
                if verbose:
                    print("first line fragment: plot from", start_horizontal_line, "to", end_horizontal_line_fragment)
                try:
                    array_hic[start_horizontal_line:end_horizontal_line_fragment, v_line_y_position] = val_max
                except IndexError:
                    if verbose:
                        print("uuh there was an error", start_horizontal_line, end_horizontal_line_fragment)

            # plot last vertical line fragment of line
            elif index_interval == n_subtraction_intervals - 1:

                end_previous_interval = list_subtraction_intervals[index_interval - 1][-1] + 1  # +1 needed otherwise there's clipping
                end_horizontal_line_fragment = interval[0]
                if verbose:
                    print("last line fragment: plot from", end_previous_interval, "to", end_horizontal_line_fragment)
                try:
                    array_hic[end_previous_interval:end_horizontal_line_fragment, v_line_y_position] = val_max
                except IndexError:
                    if verbose:
                        print("uuh there was an error", end_previous_interval, end_horizontal_line_fragment)

            # plot vertical line fragment of line if it's in between
            else:
                end_previous_interval = list_subtraction_intervals[index_interval - 1][-1] + 1  # +1 needed otherwise there's clipping
                end_horizontal_line_fragment = interval[0]
                if verbose:
                    print("line fragment #", index_interval + 1, ": plot from", end_previous_interval, "to", end_horizontal_line_fragment)
                try:
                    array_hic[end_previous_interval:end_horizontal_line_fragment, v_line_y_position] = val_max
                except IndexError:
                    if verbose:
                        print("uuh there was an error", end_previous_interval, end_horizontal_line_fragment)

    ###########################################
    # cut horizontal lines that cross a window and mark each line fragment on the array

    if verbose:
        print("\n-------------------------\nHorizontal lines")

    for key in dict_loops.keys():
        list_bin_indices = sorted(dict_loops[key], reverse=True)

        h_line_min = list_bin_indices[1]
        h_line_max = list_bin_indices[0]

        if verbose:
            print("\nVertical Line for loop #", key, "has starts and stops at x-axis bin indices:", h_line_min, h_line_max)

        h_line_y_position = h_line_max

        if verbose:
            print("x min: ", h_line_min, "x max: ", h_line_max)

        set_subtraction_intervals = dict_horizontal_line_cutouts[h_line_y_position]
        if verbose:
            print("list subrtraction intervals vor vertical line on X position =", h_line_y_position)
            print(set_subtraction_intervals)

        # start getting the horizontal lines to plot
        list_unsorted_subtraction_intervals = list(set_subtraction_intervals)
        list_subtraction_intervals = sorted(list_unsorted_subtraction_intervals, key=lambda r: r.start)
        n_subtraction_intervals = len(list_subtraction_intervals)

        for index_interval in range(len(list_subtraction_intervals)):

            interval = list_subtraction_intervals[index_interval]

            if verbose:
                print("removing inverval", interval[0], "-", interval[-1], "from vertical line at Xpos =", h_line_y_position)

            # plot first horizontal line fragment of line
            if n_subtraction_intervals == 1:

                start_horizontal_line = interval[-1] + 1
                end_horizontal_line_fragment = h_line_max
                if verbose:
                    print("first line fragment: plot from", start_horizontal_line, "to", end_horizontal_line_fragment)
                try:
                    array_hic[h_line_y_position, start_horizontal_line:end_horizontal_line_fragment] = val_max
                except IndexError:
                    if verbose:
                        print("uuh there was an error", start_horizontal_line, end_horizontal_line_fragment)

            # plot last horizontal line fragment of line
            elif index_interval == n_subtraction_intervals - 1:

                end_previous_interval = interval[-1] + 1  # +1 needed otherwise there's clipping
                end_horizontal_line_fragment = h_line_max
                if verbose:
                    print("last line fragment: plot from", end_previous_interval, "to", end_horizontal_line_fragment)
                try:
                    array_hic[h_line_y_position, end_previous_interval:end_horizontal_line_fragment] = val_max
                except IndexError:
                    if verbose:
                        print("uuh there was an error", end_previous_interval, end_horizontal_line_fragment)

            # plot horizontal line fragment of line if it's in between
            else:
                end_previous_interval = interval[-1] + 1  # +1 needed otherwise there's clipping
                end_horizontal_line_fragment = list_subtraction_intervals[index_interval + 1][0]
                if verbose:
                    print("line fragment #", index_interval + 1, ": plot from", end_previous_interval, "to", end_horizontal_line_fragment)
                try:
                    array_hic[h_line_y_position, end_previous_interval:end_horizontal_line_fragment] = val_max
                except IndexError:
                    if verbose:
                        print("uuh there was an error", end_previous_interval, end_horizontal_line_fragment)

    return array_hic


if __name__ == "__main__":
    pass
