from scipy import signal as sig
import csv
import sys
import script_generate_marked_input_vector as vecsvg
import pandas as pd
from scipy import stats
import py3dpolys_le.hic_analysis as ha
import numpy as np

"""
notes

peak calling with my methods needs very different thresholds that are not comparable at all
peak calling without threshold leads to an immense amount of peaks that's not suitable for visualization

correlation without peak-calling is very low unless very coarsly binned


"""


def read_chip_outfile(filepath_chip_outfile):

    csv.field_size_limit(sys.maxsize)  # otherwise python crashes

    print("reading outfile")

    outfile = open(filepath_chip_outfile, "r")

    readerobject_outfile = csv.reader(outfile, delimiter="\t")

    raw_signal_string = [line for line in readerobject_outfile][-1][0].strip()

    signal_separated = raw_signal_string.split("   ")

    signal_floats = []
    for value_str in signal_separated:
        try:
            signal_floats.append(float(value_str))
        except ValueError:
            pass

    return signal_floats


def read_chip_outfile_pandas(filepath_chip_outfile):

    print("read_chip_outfile_pandas")

    signal = pd.read_csv(filepath_chip_outfile, delim_whitespace=True, encoding='utf-8')
    signal = ha.chip_seq(signal, 1)

    return signal


def read_chip_bedGraph(filepath_chip_bedGraph):

    print("reading Bedgraph")

    bedgraph = open(filepath_chip_bedGraph, "r")

    readerobject_bedgraph = csv.reader(bedgraph, delimiter="\t")

    signal = [float(line[3]) for line in readerobject_bedgraph]

    return signal


def read_chip_bedGraph_pandas(filepath_chip_bedGraph):

    exp_chip_pd = pd.read_csv(filepath_chip_bedGraph, delim_whitespace=True, names=['chrom', 'start', 'end', 'value'],
                              encoding='utf-8')
    exp_chip_x_pd = exp_chip_pd[exp_chip_pd.chrom.isin(['6', 'chrX', 'X'])]
    exp_chip_out = exp_chip_x_pd[['value']].values.transpose()

    return exp_chip_out


def chip_peak_call(filepath_chip, height=0.0007):

    filetype = filepath_chip.split(sep="/")[-1].split(".")[-1]

    if filetype != "out":
        #signal = read_chip_bedGraph(filepath_chip)
        signal = read_chip_bedGraph_pandas(filepath_chip)[0]

    else:
        signal = read_chip_outfile_pandas(filepath_chip)
        #signal = read_chip_outfile(filepath_chip)

    #print(signal)
    signal_norm = signal / sum(signal)
    #peak_indices, threshold_values = sig.find_peaks(signal_norm, threshold=threshold)
    peak_indices, threshold_values = sig.find_peaks(signal_norm, height)
    #print(peak_indices, threshold_values)

    return peak_indices


def get_chip_correlation_of_simulation_to_experiment(filepath_experiment_bedGraph,
                                                     filepath_simulation_outfile,
                                                     binsize=20):

    chip_sim = pd.read_csv(filepath_simulation_outfile, delim_whitespace=True, encoding='utf-8')
    chip_sim_binned = ha.chip_seq(chip_sim, binsize)
    chip_sim_binned_normed = chip_sim_binned / sum(chip_sim_binned)

    chip_invivo = read_chip_bedGraph_pandas(filepath_experiment_bedGraph)
    chip_invivo_binned = ha.chip_seq(chip_invivo, binsize)
    chip_invivo_binned_normed = chip_invivo_binned / sum(chip_invivo_binned)

    print("sim chip - length =", len(chip_sim_binned_normed), ": ", chip_sim_binned_normed[:10])
    print("exp chip - length =", len(chip_invivo_binned_normed), ": ", chip_invivo_binned_normed[:10])

    correlation = stats.spearmanr(chip_invivo_binned_normed, chip_sim_binned_normed)

    return correlation


if __name__ == "__main__":

    filepath_chip_exp = "../../py3dpolys_le/data/ce/brejc/DPY27_N2_L3_average_ce11_2kb_chrX.bedGraph"
    filepath_chip_sim = "../../py3dpolys_le/data/simulations/chip_seqs/hybrid_l_rex_b_rex_bid_sym.chip.out"

    filepath_chromosome_sizes = "../../py3dpolys_le/data/ce/chromosome_sizes_chrIVXnotation.tsv"

    binsize = 20

    correlation = get_chip_correlation_of_simulation_to_experiment(filepath_chip_exp, filepath_chip_sim, binsize=binsize)
    print("correlation =", correlation)

    folderpath_sim = "/".join(filepath_chip_sim.split(sep="/")[:-1])+"/"
    filename_sim = filepath_chip_sim.split(sep="/")[-1].split(sep=".")[0]

    correlationfile = open(folderpath_sim + "correlation_chip_" + filename_sim + "_v_invivo." + str(binsize) + ".txt", "w")
    correlationfile.write(str(correlation[0]))
    correlationfile.close()

    height = 0.0007
    peaks_sim = chip_peak_call(filepath_chip_sim, height=height)
    peaks_sim_genomic = peaks_sim * 2000
    print(len(peaks_sim_genomic))
    vecsvg.draw_input_vector_svg(peaks_sim_genomic, filepath_chromosome_sizes, folderpath_sim + "chip_" + filename_sim + "_height-" + str(height))

    peaks_exp = chip_peak_call(filepath_chip_exp, height=height)
    peaks_exp_genomic = peaks_exp * 2000
    print(len(peaks_exp_genomic))
    vecsvg.draw_input_vector_svg(peaks_exp_genomic, filepath_chromosome_sizes, "chip_invivo.height-" + str(height))

    """
    import boundary_files as bdf
    import math
    length = math.ceil(bdf.read_chromosome_sizes(filepath_chromosome_sizes)["chrX"] / 2000)
    vector_peaks_sim = np.zeros(length)
    vector_peaks_sim[peaks_sim] = 1
    print("vetor peaks sim", vector_peaks_sim)

    vector_peaks_exp = np.zeros(length)
    vector_peaks_exp[peaks_exp] = 1
    print("vetor peaks sim", vector_peaks_exp)

    correlation = stats.spearmanr(vector_peaks_exp, vector_peaks_sim)

    print(correlation)
    """
