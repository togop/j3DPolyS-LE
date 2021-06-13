import csv


def bed_to_loadingfile(filepath_input_bed, filepath_output_loadingfile):

    boundaryfile = open(filepath_input_bed, "r")

    readerobject_boundaryfile = csv.reader(boundaryfile, delimiter="\t")

    boundaryfile_output = open(filepath_output_loadingfile, "w")

    boundaryfile_output.write("name,position,length,probability\n")

    for line in readerobject_boundaryfile:
        print(line)
        name = line[0] + "_" + line[1]
        midpoint = (int(line[1]) + int(line[2])) / 2

        line_output = name + "," + str(round(midpoint)) + ",1,1\n"

        print(line_output)
        boundaryfile_output.write(line_output)


def run_test_bed_to_loadingfile():

    filepath_input_bed = "/media/cubix/D86E-6C50/loops/published/N2.chrX.allValidPairs.hic.5-10kbAnchors.bed"
    filepath_output_loadingfile = "/media/cubix/D86E-6C50/boundaries/sip-loop-anchors-loading.csv"

    bed_to_loadingfile(filepath_input_bed,filepath_output_loadingfile)


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


if __name__ == "__main__":
    run_test_bed_to_loadingfile()