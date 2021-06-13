import csv


def read(filepath_chromosome_sizes):
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


if __name__ == "__main__":
    pass
