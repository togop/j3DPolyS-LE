#! /usr/bin/env python

import argparse
import csv
import logging
import os
import re

import sys
from filelock import SoftFileLock

from . import hic_analysis as ha
from . import plot_hic
# from _version import __name__

dummy_sim = False  # set to False; True = dummy simulation mode: echo commands only


def main():

    # Initialization
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)
    logger = logging.getLogger('3dpolys_le_stats')

    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output_folder", default=".",
                   help=f"Simulation's output folder containing "
                        f"result files: {ha.CONFIG_OUT}, {ha.CONTACT_OUT}, {ha.DR_OUT}.")
    p.add_argument("-a", "--analyse", default="./out/analyse", help="'Analyse' step output folder containing hic_*.hdf5 files.")
    p.add_argument("-b", "--boundary", default=None, help="Boundary file in TSV format.")
    p.add_argument("-bd", "--boundary_direction", default=0, help="Impermeability direction applied to all boundaries: "
                                                                  "-1: opposite direction, 0: both, 1: same direction.")
    p.add_argument("-bf", "--boundary_factor", default=1., help="Impermeability factor applied to all boundaries.")
    p.add_argument("-bs", "--boundary_score", help="Whether individual boundary scores were applied in a simulation.",
                   action='store_true')
    p.add_argument("-t", "--tads_boundary", default=None,
                   help="TADs boundary file in CSV format (same as boundary.csv) to be used for calculating chi2-min score. "
                        "Also supported Loops file in .bed.tsv format, with the following columns: "
                        "chromosome  anchor1  anchor2")
    p.add_argument("-e", "--exp_cool", default=f"{ha.PUBLISHED_FOLDER}/N2_5000b.cool",
                   help="Experimental cool file with which simulation data to be compared.")
    p.add_argument("-l", "--nlef", default="200", help="Nlef value used in a simulation.")  # TODO optional read it form input.dat
    p.add_argument("-m", "--km", default="2.7e-3", help="km value used in a simulation.")  # TODO optional read it form input.dat
    p.add_argument("-i", "--input_cfg", default="./input.cfg", help="input.cfg file used in a simulation.")
    p.add_argument("-r", "--radius_contact", default=0., help="Contact radius in lattice units (1=70nm) "
                                                              "used in the 'analyse' step to extract Hi-C matrixes.",
                   type=float)
    p.add_argument("-cp", "--contact_probability", help="In combination with the contact radius parameter, "
                                                        "whether to use contact radius probability with the formula:"
                                                        " (1 - r^2 / max_r^2), "
                                                        "where <max_r> is the value of the --radius_contact parameter."
                   , action='store_true')
    p.add_argument("-f", "--stats_file", default="./sim_stats.csv", help="Simulation statistics' repository file.")
    p.add_argument("--bin_size", default=1, help="Chip-seq bin size used for plotting. Default: 1 = 2kb.", type=int)
    p.add_argument("--replace", help="Whether to replace existing files.", action='store_true')
    p.add_argument("--chi2_mode", default=ha.CHI2_MODE_LINEAR, choices=['log', 'linear'],
                   help="Chi2-min mode for sampling contact distances.")
    p.add_argument("--cmp_chrs", nargs='*', default=None,
                   help="Synonyms of the Hi-C chromosome with which simulation data to be compared. "
                        "Default: hic_analysis.py::CHR_SYNONYMS!")

    args = p.parse_args(sys.argv[1:])

    logger.info(f'start with parameters: {args}')

    if args.cmp_chrs is not None:
        logger.warning(f'Overwriting the default hic_analysis.CHR_SYNONYMS: { ",".join(args.cmp_chrs)} '
                       f'with which HiCs will be compared!')
        ha.CHR_SYNONYMS = args.cmp_chrs

    resolutions = [int(i * ha.EXP_RESOLUTION) for i in ha.EXP_FACTORS]

    exp_cool = args.exp_cool
    boundary = args.boundary
    tads_boundary = args.tads_boundary

    cmp_hic_file = re.sub(r'.cool', '', exp_cool)  # hic to compare with

    # find last HIC.hdf5 do analysis and store
    sim_hic_file = ha.get_last_hic(args.analyse)

    # need only normed for chi2_log and chi2_linear and for given tads-boundary sites and 1tad(the whole chromosome)
    (chi2_lin, alpha_lin) = ha.compare_hic_chromosome(sim_hic_file, cmp_hic_file, chrs=ha.CHR_SYNONYMS, res=ha.RESOLUTION,
                                                      tads_boundary=tads_boundary, norm=True,
                                                      chi2_mode=ha.CHI2_MODE_LINEAR)  # , plot=True)
    (chi2_log, alpha_log) = ha.compare_hic_chromosome(sim_hic_file, cmp_hic_file, chrs=ha.CHR_SYNONYMS, res=ha.RESOLUTION,
                                                      tads_boundary=tads_boundary, norm=True, chi2_mode=ha.CHI2_MODE_LOG)

    # import matplotlib.pyplot as plt
    # plt.plot(exp_score, sim_score, '.')  # juts to visualize with what the correlation coefficient has to deal with

    logger.info(
        f'comparing {sim_hic_file} with {exp_cool} result chi2: {chi2_log}')

    # save stats in repository file
    if not os.path.isfile(args.stats_file):
        # header
        with open(args.stats_file, mode='w') as stats_file:
            stats_writer = csv.writer(stats_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            # !!! it is important the first column to be unique as this will be used by panda as index !!!
            stats_writer.writerow(
                ['sim_hic_file', 'sim_out_folder', 'exp_cool', 'resolution',
                 'boundary', 'boundary_direction', 'boundary_factor', 'boundary_score',
                 'tads_boundary', 'input.cfg', 'nlef', 'km', 'radius_contact', 'chi2_log', 'alpha_log',
                 'chi2_lin', 'alpha_lin', 'chr'])

    with SoftFileLock(f'{args.stats_file}.lock'):
        with open(args.stats_file, mode='a+') as stats_file:
            stats_writer = csv.writer(stats_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            cp = 'p' if args.contact_probability else ''
            radius_p = f'{args.radius_contact}{cp}' if args.radius_contact > 0 else ''
            stats_writer.writerow(
                [sim_hic_file, args.output_folder, exp_cool, ha.RESOLUTION,
                 boundary, args.boundary_direction, args.boundary_factor, args.boundary_score,
                 tads_boundary, args.input_cfg, args.nlef, args.km, radius_p, chi2_log, alpha_log,
                 chi2_lin, alpha_lin, ha.CHR_SYNONYMS[-1]])

    # generate hic_*_hot_r.png plot files if not already done
    if not os.path.exists(re.sub('.hdf5', '_hot_r.png', sim_hic_file)):
        plot_hic.run(args.analyse, cmap="hot_r", file_ext="png")


if __name__ == '__main__':
    main()
