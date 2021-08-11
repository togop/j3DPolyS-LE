import h5py
import cooler
import numpy as np
import pandas as pd
from PIL import Image
from PIL import ImageEnhance
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter

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


def info_old(array_hic, verbose=False):
    # get max/min value in method array
    val_max = np.max(array_hic)
    val_min = np.min(array_hic)
    shape = array_hic.shape
    if verbose:
        print("Hi_C:: max", val_max, " / min", val_min, "| Shape: ", shape)

    return val_max, val_min


def info(array_hic, arrayname=None):
    val_max = np.max(array_hic)
    val_min = np.min(array_hic)
    shape = array_hic.shape
    print("Hi_C:", str(arrayname), ": max", val_max, " / min", val_min, "| Shape: ", shape)


def get_max_min(array_hic):
    val_max = np.max(array_hic)
    val_min = np.min(array_hic)
    return val_max, val_min


def plot_new(hic_array, filename_output=None, set_clim=None, save=False, resolution=2000):

    hic_array = hic_array.astype("float")

    # make image window
    fig = plt.figure()

    plt.imshow(hic_array, interpolation='nearest', cmap="hot_r")
    if set_clim:
        plt.clim(set_clim)

    plt.colorbar()

    plt.gca().get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: int(x * resolution)))
    plt.gca().get_yaxis().set_major_formatter(FuncFormatter(lambda y, p: int(y * resolution)))

    if filename_output is not None:
        plt.title(filename_output)
        if save:
            fig.savefig(filename_output, dpi=3000)


def plot_old(hic_array, filename_output=None, set_clim=(-2.75, 0)):

    val_max, val_min = info_old(hic_array)

    # show hic image
    fig = plt.figure()

    plt.imshow(hic_array, interpolation='nearest', cmap="Greys")
    plt.clim(set_clim)

    plt.colorbar()
    plt.title(filename_output)

    if filename_output is not None:
        fig.savefig(filename_output, dpi=2000)


def plot_sim_for_figures(hic_array, filename_output=None, set_clim=(-2.75, 0), resolution=2000):

    hic_array = hic_array.astype("float")

    val_max, val_min = info_old(hic_array)

    # show hic image
    fig = plt.figure()

    plt.imshow(hic_array, interpolation='nearest', cmap="Greys")
    plt.clim(set_clim)

    lambda_x = lambda x, p: str((x * resolution) / 1000000).split(sep=".")[0] + "mb"
    lambda_y = lambda y, p: str((y * resolution) / 1000000).split(sep=".")[0] + "mb"
    plt.gca().get_xaxis().set_major_formatter(FuncFormatter(lambda_x))
    plt.gca().get_yaxis().set_major_formatter(FuncFormatter(lambda_y))

    if filename_output is not None:
        fig.savefig(filename_output, dpi=1000)


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


def remove_nan_from_array(hic_array):

    indicex_nan = (np.argwhere(np.isnan(hic_array)))

    for index in indicex_nan:
        y = index[0]
        x = index[1]
        hic_array[y, x] = 0

    return hic_array


def remove_ninf_from_array(input_array):

    indicex_nan = np.argwhere(input_array == np.NINF)

    # this is super slow
    # hic_array[indicex_nan] = 0

    # this is fast but needs 8.33 TiB of ram
    # shape = input_array.shape
    # output_array = np.zeros((shape))
    # output_array[indicex_nan] = input_array[indicex_nan]

    # this is slow, but mindnumbingly so
    for index in indicex_nan:
        y = index[0]
        x = index[1]
        input_array[y, x] = 0

    return input_array


def change_contrast_IDONTTHINKTHISWORKS(hic_array, factor, verbose=False):
    minval = np.percentile(hic_array, 2)
    maxval = np.percentile(hic_array, 99.999)

    max_ = np.max(hic_array)
    min_ = np.min(hic_array)

    hic_array = np.clip(hic_array, minval, maxval)
    max_new = np.max(hic_array)
    min_new = np.min(hic_array)
    hic_array = ((hic_array - minval) / (hic_array - minval))

    # hic_array = ((hic_array - hic_array.min()) / (hic_array.max() - hic_array.min())) * 1

    return hic_array


def contrast(array_image, contrast_adjustment, bitrange=256):

    # get array data
    arrshape = array_image.shape
    arrmin = np.min(array_image)
    arrmax = np.max(array_image)

    # get bitrange data
    image_max = range(bitrange)[-1]
    image_half = int(bitrange / 2)
    image_min = range(bitrange)[0]

    # calculate contrast factor
    a = image_max * (contrast_adjustment + image_max)
    b = image_max * (image_max - contrast_adjustment)
    contrast_factor = a / b

    # transform array
    for y in range(arrshape[0]):
        for x in range(arrshape[1]):
            originalvalue = array_image[y, x]
            newvalue = (contrast_factor * (originalvalue - image_half)) + image_half
            array_image[y, x] = newvalue

    # clip array to given bitrange
    array_image = np.clip(array_image, image_min, image_max)

    return array_image


def stretch_value_range(array_input, bitrange=256):

    # get data from array
    array_input = array_input.astype(np.uint8)  # the array should actually be integer to begin with
    maxval_input = round(array_input.max())
    minval_input = round(array_input.min())

    length_value_range_input = round((maxval_input - minval_input) + 1)
    maxval_input_plusone = maxval_input + 1

    # get data for given bitrange
    lowest_bitrange_value = range(bitrange)[0]
    highest_bitrange_value = range(bitrange)[-1]

    # generate look up table
    LUT = np.zeros(bitrange, dtype=np.uint8)
    LUT[minval_input:maxval_input_plusone] = np.linspace(start=lowest_bitrange_value, stop=highest_bitrange_value, num=length_value_range_input, endpoint=True, dtype=np.uint8)

    array_output = LUT[array_input]

    return array_output


def transform_array_newmax(array_input, int_newmax, custom_oldmax=None):

    if custom_oldmax:
        factor_for_multiplication = float(int_newmax / custom_oldmax)

    else:
        max_inputarray = np.max(array_input)
        factor_for_multiplication = float(int_newmax / max_inputarray)

    output_array = np.multiply(array_input, factor_for_multiplication)
    if custom_oldmax:
        output_array = np.clip(output_array, 0, int_newmax)

    return output_array


def change_contrast_PILLOW(hic_array, factor, verbose=False):
    if verbose:
        val_max, val_min = info_old(hic_array, True)
    else:
        val_max, val_min = info_old(hic_array)

    # make numbers in array into a usable range for pillow (l image mode has pixels from 0 to 255)
    factor_for_multiplication = 255 / val_max
    hic_array = np.multiply(hic_array, factor_for_multiplication)

    hic_image = Image.fromarray(hic_array)
    hic_image = hic_image.convert("L")

    enhancer = ImageEnhance.Contrast(hic_image)

    hic_enhanced = enhancer.enhance(factor)

    hic_array_output = np.array(hic_enhanced)

    return hic_array_output


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


def normalize_invivo_array(hic_array, alpha=1.0, verbose=False):



    if verbose:
        print("\nstart normalizing invivo hic array with proprieties:")
        unused_var_1, unused_var_2 = info_old(hic_array, True)
        print("remove diagonal from hic array")

    hic_array = hic_array * (255 / unused_var_1)

    # fill the empty middle lane with max value in array
    #np.fill_diagonal(hic_array, unused_var_1)

    # get data from array
    if verbose:
        val_max, val_min = info_old(hic_array, True)
    else:
        val_max, val_min = info_old(hic_array)

    # divide matrix so it has val_max = 1
    if verbose:
        print("divide array by factor", val_max)
    hic_array = np.divide(hic_array, val_max)
    #histogram(hic_array)
    #hic_array = np.divide(hic_array, alpha * 10)

    # log transform matrix
    if verbose:
        unused_var_1, unused_var_2 = info_old(hic_array, True)
        print("apply log10 transform on array")
    hic_array = np.log10(hic_array)

    if verbose:
        print("\ndone normalizing invivo hic array")
        unused_var_1, unused_var_2 = info_old(hic_array, True)
        print("-------------------")

    return hic_array


def normalize_simulation_array(hic_array, verbose=False):

    if verbose:
        print("\nstart normalizing simulation hic array with proprieties:")
        unused_var_1, unused_var_2 = info_old(hic_array, True)
        print("remove diagonal from hic array")

    # fill the empty middle lane with max value in array
    np.fill_diagonal(hic_array, 0)

    # get data from array
    if verbose:
        val_max, val_min = info_old(hic_array, True)
    else:
        val_max, val_min = info_old(hic_array)

    # divide matrix so it has val_max = 1
    if verbose:
        print("divide array by factor", val_max)
    hic_array = np.divide(hic_array, val_max)

    # log transform matrix
    if verbose:
        unused_var_1, unused_var_2 = info_old(hic_array, True)
        print("apply log10 transform on array")
    hic_array = np.log10(hic_array)

    if verbose:
        print("\ndone normalizing simulation hic array")
        unused_var_1, unused_var_2 = info_old(hic_array, True)
        print("-------------------")

    return hic_array


def get_diagonal_as_list(array, offset=1):

    y_index = offset
    x_index = 0

    maxindex = array.shape[0] - 1
    list = []
    while y_index <= maxindex:
        list.append(array[y_index, x_index])
        y_index += 1
        x_index += 1

    return list


def get_median_of_diagonal(array_input, diagonal_thickness, verbose=False):

    diagonaloffset_start = 1
    diagonaloffset_end = diagonal_thickness + 1

    if verbose:
        print(range(diagonaloffset_start, diagonaloffset_end))

    list_diagonal = []
    for offset in range(diagonaloffset_start, diagonaloffset_end):

        if verbose:
            print(offset)

        list_diagonal = list_diagonal + get_diagonal_as_list(array_input, offset)

    median = np.median(list_diagonal)

    return median


def merge_two_hic_arrays(hic_array_upper, hic_array_lower, verbose=False):

    # get data from array
    if verbose:
        val_max_upper, val_min_upper = info_old(hic_array_upper, True)
        val_max_lower, val_min_lower = info_old(hic_array_lower, True)
        val_max = max([val_max_upper, val_max_lower])
    else:
        val_max_upper, val_min_upper = info_old(hic_array_upper)
        val_max_lower, val_min_lower = info_old(hic_array_lower)
        val_max = max([val_max_upper, val_max_lower])

    # delete bottom half of hic-array
    for index_row in range(len(hic_array_upper)):
        hic_array_upper[index_row][:index_row] = hic_array_lower[index_row][:index_row]

    # fill the empty middle lane with max value in array
    np.fill_diagonal(hic_array_upper, 0)

    return hic_array_upper


def resize(hic_array, shape_array):

    hic_array = np.multiply(hic_array, 1000)
    image_hic = Image.fromarray(hic_array)
    image_hic = image_hic.convert("L")
    image_resized = image_hic.resize(shape_array)
    hic_array = np.array(image_resized)

    return hic_array


def histogram(normalized_hic_array):

    normalized_hic_array[np.where(normalized_hic_array == np.NINF)] = -2.75
    normalized_hic_array = np.add(normalized_hic_array, 2.75)
    normalized_hic_array = np.multiply(normalized_hic_array, 100)
    flat = normalized_hic_array.flatten()
    plt.hist(flat, bins=list(range(0, 275, 2)))


def histogram_new(array_input, verbose=False, filename=None):

    arrmax, arrmin = get_max_min(array_input)
    histogram_bins = np.linspace(arrmin, arrmax, 100)

    if verbose:
        print("generating histogram of array with max", arrmax, "and min", arrmin)

    flattened_array = array_input.flatten().tolist()
    if verbose:
        print("array flatened")

    flattened_array = [i for i in flattened_array if i != arrmin]
    if verbose:
        print("filtering done")

    # make hist plot
    fig = plt.figure()
    plt.hist(flattened_array, bins=histogram_bins)
    if filename:
        fig.savefig(filename)


def save_with_pillow(hic_array):

    image_hic = Image.fromarray(hic_array)
    image_hic.save("hic_array.tiff")


def mark_loop_windows(array_hic, filepath_boundaryfile_csv, resolution_loops, verbose=False):

    if verbose:
        print("\nstart marking loop windows on hic array")

    # get data from array
    if verbose:
        print("input array data:")
        val_max, val_min = info_old(array_hic, True)
        shape_array_hic = np.shape(array_hic)
        index_max_array_hic = shape_array_hic[0]

    else:
        val_max, val_min = info_old(array_hic, False)
        shape_array_hic = np.shape(array_hic)
        index_max_array_hic = shape_array_hic[0]

    # mark anchors and loops
    thickness_line = 1
    wideness_window = 6
    trace_color = val_max  # 2.4065402  # 0

    # get loop coordinates
    dict_loops = loops.get_dict_from_boundaryfile(filepath_boundaryfile_csv, resolution=resolution_loops)
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
        array_hic[window_min_y, window_min_x:window_max_x] = trace_color

        # plot bot window line
        array_hic[window_max_y, window_min_x:window_max_x + 1] = trace_color  # the 1 here is needed for nice squares

        # plot left window line
        array_hic[window_min_y:window_max_y, window_min_x] = trace_color

        # plot left window line
        array_hic[window_min_y:window_max_y, window_max_x] = trace_color

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
                    array_hic[start_horizontal_line:end_horizontal_line_fragment, v_line_y_position] = trace_color
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
                    array_hic[end_previous_interval:end_horizontal_line_fragment, v_line_y_position] = trace_color
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
                    array_hic[end_previous_interval:end_horizontal_line_fragment, v_line_y_position] = trace_color
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
                    array_hic[h_line_y_position, start_horizontal_line:end_horizontal_line_fragment] = trace_color
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
                    array_hic[h_line_y_position, end_previous_interval:end_horizontal_line_fragment] = trace_color
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
                    array_hic[h_line_y_position, end_previous_interval:end_horizontal_line_fragment] = trace_color
                except IndexError:
                    if verbose:
                        print("uuh there was an error", end_previous_interval, end_horizontal_line_fragment)

    return array_hic


def find_lowest_nonzero_value(array_input):

    percentile_sequence = np.arange(1, 101, 1)

    for percentile in percentile_sequence:
        value = np.percentile(array_input, percentile)
        #print(value)
        if value > 0:
            return value


if __name__ == "__main__":
    pass


