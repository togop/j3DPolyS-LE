import csv
import math
import pandas as pd


def read(filepath_input_boundaryfile, bool_return_header=False):
    boundaryfile = open(filepath_input_boundaryfile, "r")

    # create boundary file array
    readerobject_boundaryfile = csv.reader(boundaryfile, delimiter=",")
    array_boundaries = [line[:3] for line in readerobject_boundaryfile]  # [:3] so this only reads the first 3 columns
    header_boundaryfile = ','.join(array_boundaries[0])  # get header
    array_boundaries = array_boundaries[1:]  # crop header
    boundaryfile.close()

    if bool_return_header:
        return array_boundaries, header_boundaryfile

    else:
        return array_boundaries


def read_chromosome_sizes(filepath_chromosome_sizes):
    chromosome_sizes_file = open(filepath_chromosome_sizes, "r")

    # create chromosome sizes array
    readerobject_chromosome_sizes = csv.reader(chromosome_sizes_file, delimiter="\t")
    chromosome_sizes = [line for line in readerobject_chromosome_sizes]
    chromosome_sizes_file.close()

    dict_chromosome_sizes = {}
    for chromosome in chromosome_sizes:
        chr_name = chromosome[0]
        chr_length = int(chromosome[1])
        dict_chromosome_sizes[chr_name] = chr_length

    return dict_chromosome_sizes


def boundaryfile_to_loadingfile(fielpath_input_boundaryfile, filepath_output_loadingfile):

    boundaryfile = open(fielpath_input_boundaryfile, "r")

    readerobject_boundaryfile = csv.reader(boundaryfile, delimiter=",")

    boundaryfile_output = open(filepath_output_loadingfile, "w")

    boundaryfile_output.write("name,position,length,probability\n")

    readerobject_boundaryfile.__next__()  # ignore header

    for line in readerobject_boundaryfile:
        print(line)
        name = line[0]
        midpoint = int(line[1])

        line_output = name + "," + str(midpoint) + ",1,1\n"

        print(line_output)
        boundaryfile_output.write(line_output)


def run_test_boundaryfile_to_loadingfile():

    filepath_boundary_sites = "/media/cubix/D86E-6C50/boundaries/3dpolysle_data/boundary_sites/boundaries_chrX_5-10kbLoops_chromosight_auto.csv"

    filepath_loading_sites = "/media/cubix/D86E-6C50/boundaries/3dpolysle_data/loading_sites/loading_sites_chromosight.csv"

    boundaryfile_to_loadingfile(filepath_boundary_sites, filepath_loading_sites)


def adjust_boundaries_for_X_fusion_chromatin(array_boundaryfile_rows, length_of_chromosome_fused_to_X, resolution_boundaryfile=5000, verbose=False):

    value_to_add_to_boundaries = math.ceil(length_of_chromosome_fused_to_X / resolution_boundaryfile) * resolution_boundaryfile

    if verbose:
        print(value_to_add_to_boundaries)

    array_adjusted_boundaryfile_rows = []
    for row in array_boundaryfile_rows:

        if verbose:
            print("\n")
            print(row)

        adjusted_boundarysite = int(row[1]) + value_to_add_to_boundaries
        adjusted_bpos = round(adjusted_boundarysite / resolution_boundaryfile)

        if verbose:
            print(adjusted_boundarysite - value_to_add_to_boundaries)
            print(adjusted_boundarysite - int(row[1]))

        line = [row[0], str(adjusted_boundarysite), row[2], row[3], str(adjusted_bpos), row[5]]

        if verbose:
            print(line)

        array_adjusted_boundaryfile_rows.append(line)

    return array_adjusted_boundaryfile_rows


def write(array_boundaryfile_rows, filepath_output_boundaryfile):

    output_boundaryfile = open(filepath_output_boundaryfile, "w")

    output_boundaryfile.write('name,midpoint,impermeability,score,b-position,strand\n')

    for row in array_boundaryfile_rows:

        output_line = row[0] + "," + row[1] + "," + row[2] + "," + row[3] + "," + row[4] + "," + row[5] + "\n"
        output_boundaryfile.write(output_line)


def subtract_thisdoesnotworkatall(boundaries_base, boundaries_to_subtract, resolution=2000):

    for index_subtraction_boundary in range(len(boundaries_to_subtract)):
        for index_base_boundary in boundaries_base:
            bpos_base = round(int(index_base_boundary[1]) / resolution)
            bpos_subtraction = round(int(index_subtraction_boundary[1]) / resolution)
            if bpos_subtraction == bpos_base:
                print(bpos_subtraction, bpos_base)
                print("boundary", index_subtraction_boundary, "was subtracted from base", index_base_boundary)
                boundaries_base.remove(index_subtraction_boundary)


if __name__ == "__main__":

    """
    # test boundary subtraction
    filepath_base_boundaries = "../data/ce/boundary_sites/boundaries_chrX_manual_detection.csv"
    filepath_subtraction_boundaries = "../data/ce/boundary_sites/boundaries_rex_sites.csv"

    boundaries_base = read(filepath_base_boundaries)
    boundaries_to_subtract = read(filepath_subtraction_boundaries)

    subtract_thisdoesnotworkatall(boundaries_base, boundaries_to_subtract)
    """

    pass


