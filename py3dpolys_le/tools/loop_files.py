# python channel libraries
import csv


def get_dict_from_boundaryfile(filepath_csv, resolution=None, header=True):

    file_csv = open(filepath_csv, "r")

    readerobject_csv = csv.reader(file_csv, delimiter=",")

    if header:
        next(readerobject_csv)  # skip header

    dict_loops = {}
    index_loop = 0
    anchor_counter = 1
    list_loop_anchors = []
    for line in readerobject_csv:

        anchor_basepair_position = int(line[1])

        if resolution is not None:  # return bin position instead if resolution is provided
            anchor_basepair_position = round(anchor_basepair_position / resolution)

        list_loop_anchors.append(anchor_basepair_position)

        if anchor_counter >= 2:
            dict_loops[index_loop] = list_loop_anchors
            index_loop += 1
            anchor_counter = 1
            list_loop_anchors = []
        else:
            anchor_counter += 1

    file_csv.close()

    return dict_loops


def loops_dict_from_bedpe(filepath_bedpe, resolution=None, header=False):

    file_bedpe = open(filepath_bedpe, "r")

    readerobject_bedpe = csv.reader(file_bedpe, delimiter="\t")

    if header:
        next(readerobject_bedpe)  # skip header

    dict_loops = {}

    index_loop = 0
    for line in readerobject_bedpe:

        # read anchor 1 data
        anchor_1_start = int(line[1])
        anchor_1_end = int(line[2])
        anchor_1_midpoint = round((anchor_1_start + anchor_1_end) / 2)

        # read anchor 2 data
        anchor_2_start = int(line[4])
        anchor_2_end = int(line[5])
        anchor_2_midpoint = round((anchor_2_start + anchor_2_end) / 2)

        if resolution is not None:  # return bin position instead if resolution is provided
            anchor_1_midpoint = round(anchor_1_midpoint / resolution)
            anchor_2_midpoint = round(anchor_2_midpoint / resolution)

        # write dict
        dict_loops[index_loop] = [anchor_1_midpoint, anchor_2_midpoint]

        # go to next loop
        index_loop += 1

    return dict_loops


def run_test_loop_readers():
    filepath_manual_detection_csv = "/media/cubix/D86E-6C50/boundaries/manual_loop_anchors_boundaries.csv"
    filepath_sip_csv = "/media/cubix/D86E-6C50/boundaries/sip_loopanchor_boundaries.csv"
    filepath_bedpe = "/media/cubix/D86E-6C50/loops/sip/N2.chrX.allValidPairs.hic.5-10kbLoops.bedpe"
    dict_csv_method = get_dict_from_boundaryfile(filepath_sip_csv)  # , resolution=2000)
    dict_bedpe_method = loops_dict_from_bedpe(filepath_bedpe)  # , resolution=2000)

    print(dict_csv_method)
    print(dict_bedpe_method)


if __name__ == "__main__":
    run_test_loop_readers()
