#! /usr/bin/env python

import argparse
import csv
import logging
import os
import re

import sys
from filelock import SoftFileLock

from py3dpolys_le import hic_analysis as ha
from py3dpolys_le import plot_hic
from py3dpolys_le.job_runner import CfgJobRunner, convert_dat_to_cfg

# from _version import __name__

dummy_sim = False  # set to False; True = dummy simulation mode: echo commands only


def cli_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output_folder", default=".",
                   help=f"Simulation's output folder containing "
                        f"result files: {ha.CONFIG_OUT}, {ha.CONTACT_OUT}, {ha.DR_OUT}.")
    p.add_argument("-a", "--analyse", default="./out/analyse", help="'Analyse' step output folder containing hic_*.hdf5 files.")
    p.add_argument("-b", "--boundary", default=None, help="Boundary file in TSV format.")
    p.add_argument("-bd", "--boundary_direction", default=0, help="Impermeability direction applied to all boundaries: "
                                                                  "-1: opposite direction, 0: both, 1: same direction.")
    p.add_argument("-t", "--tads_boundary", default=None,
                   help="TADs boundary file in CSV format (same as boundary.csv) to be used for calculating chi2-min score. "
                        "Also supported Loops file in .bed.tsv format, with the following columns: "
                        "chromosome  anchor1  anchor2")
    p.add_argument("-e", "--exp_cool", default="",
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
    p.add_argument("-res", "--resolution", default=ha.SIM_RESOLUTION, type=int,
                   help=f"Resolution to downscale Chip-seq output data in the .bedGraph format. "
                        f"Default: {ha.SIM_RESOLUTION} = 2kb.")
    p.add_argument("--replace", help="Whether to replace existing files.", action='store_true')
    p.add_argument("--chi2_mode", default=ha.CHI2_MODE_LINEAR, choices=['log', 'linear'],
                   help="Chi2-min mode for sampling contact distances.")
    p.add_argument("--cmp_chrs", nargs='*', default=None,
                   help="Synonyms of the Hi-C chromosome with which simulation data to be compared. "
                        "Default: hic_analysis.py::CHR_SYNONYMS!")
    return p


def main():

    # Initialization
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)
    logger = logging.getLogger('3dpolys_le_stats')

    args = cli_parser().parse_args(sys.argv[1:])

    input_cfg = args.input_cfg
    if input_cfg.endswith('.dat'):
        input_cfg = convert_dat_to_cfg(input_cfg)

    cfg_job_runner = CfgJobRunner(input_cfg=input_cfg)

    logger.info(f'start with parameters: {args}')

    if args.cmp_chrs is not None:
        logger.info(f'Overwriting the default hic_analysis.CHR_SYNONYMS: { ",".join(args.cmp_chrs)} '
                       f'with which HiCs will be compared!')
        ha.CHR_SYNONYMS = args.cmp_chrs
    else:
        ha.CHR_SYNONYMS = re.split('\\s*;\\s*|\\s*,\\s*|\\s+', cfg_job_runner.get_sim_property(name='cmp_chrs'))

    exp_cool = args.exp_cool if args.exp_cool else cfg_job_runner.get_sim_property(name='exp_cool')
    interaction_sites = cfg_job_runner.get_sim_property(name='interaction_sites')
    lef_loading_sites = cfg_job_runner.get_sim_property(name='lef_loading_sites')
    basal_loading_factor = cfg_job_runner.get_sim_property(name='basal_loading_factor')
    boundary = args.boundary if args.boundary else cfg_job_runner.get_sim_property(name='boundary')
    tads_boundary = args.tads_boundary if args.tads_boundary else cfg_job_runner.get_sim_property(name='tads_boundary')

    cmp_hic_file = re.sub(r'\.cool|\.mcool', '', exp_cool) if exp_cool else None  # hic to compare with

    # find last HIC.hdf5 do analysis and store
    sim_hic_file = ha.get_last_hic(args.analyse)

    hic_folder = os.path.dirname(sim_hic_file)
    plots_folder = os.path.join(hic_folder, ha.PLOTS_FOLDER)

    plot_cmap = cfg_job_runner.get_property(profile='stats', name='plot_cmap', default=ha.CMAP)
    plot_format = cfg_job_runner.get_property(profile='stats', name='plot_format', default=ha.PLOT_FORMAT)

    # need only normed for chi2_log and chi2_linear and for given tads-boundary sites and 1tad(the whole chromosome)
    if exp_cool:
        (chi2_lin, alpha_lin) = ha.compare_hic_chromosome(sim_hic_file, cmp_hic_file, chrs=ha.CHR_SYNONYMS,
                                                          res=ha.RESOLUTION, tads_boundary=tads_boundary, norm=True,
                                                          chi2_mode=ha.CHI2_MODE_LINEAR,
                                                          plot_cmap=plot_cmap, plot_format=plot_format,
                                                          interaction_sites=interaction_sites,
                                                          lef_loading_sites=lef_loading_sites,
                                                          basal_loading_factor=basal_loading_factor,
                                                          lef_boundaries=boundary)
        (chi2_log, alpha_log) = ha.compare_hic_chromosome(sim_hic_file, cmp_hic_file, chrs=ha.CHR_SYNONYMS,
                                                          res=ha.RESOLUTION, tads_boundary=tads_boundary,
                                                          plots_folder=plots_folder, norm=True,
                                                          chi2_mode=ha.CHI2_MODE_LOG,
                                                          plot_cmap=plot_cmap, plot_format=plot_format,
                                                          interaction_sites=interaction_sites,
                                                          lef_loading_sites=lef_loading_sites,
                                                          basal_loading_factor=basal_loading_factor,
                                                          lef_boundaries=boundary)
    else:
        (chi2_lin, alpha_lin) = (0, 1)
        (chi2_log, alpha_log) = (0, 1)
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
                 'boundary', 'boundary_direction',
                 'tads_boundary', 'input.cfg', 'nlef', 'km', 'radius_contact', 'chi2_log', 'alpha_log',
                 'chi2_lin', 'alpha_lin', 'chr'])

    with SoftFileLock(f'{args.stats_file}.lock'):
        with open(args.stats_file, mode='a+') as stats_file:
            stats_writer = csv.writer(stats_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            cp = 'p' if args.contact_probability else ''
            radius_p = f'{args.radius_contact}{cp}' if args.radius_contact > 0 else ''
            stats_writer.writerow(
                [sim_hic_file, args.output_folder, exp_cool, ha.RESOLUTION,
                 boundary, args.boundary_direction, tads_boundary, args.input_cfg, args.nlef, args.km, radius_p,
                 chi2_log, alpha_log, chi2_lin, alpha_lin, ha.CHR_SYNONYMS[-1]])

    # generate hic_*_hot_r.png plot files if not already done
    hic_plot_file = re.sub('.hdf5', f'_{plot_cmap}.{plot_format}', sim_hic_file)
    if not os.path.exists(hic_plot_file) or args.replace:
        plot_hic.run(args.analyse, resolution=ha.SIM_RESOLUTION, cmap=plot_cmap, plot_format=plot_format)
    else:
        logger.info(f'Hic plot file already created {hic_plot_file} so skip it')

    # ChIP-seq in bedGraph format
    # chip_out_file = os.path.join(args.analyse, ha.CHIP_OUT)
    # bed_graph_file = os.path.join(args.analyse, ha.CHIP_BED_GRAPH)
    # if os.path.exists(chip_out_file) and (not os.path.exists(bed_graph_file) or args.replace):
    #    ha.chip_out_to_bedgraph(chip_out_file, bed_graph_file=bed_graph_file,
    #                            chrom=ha.SIM_CHR, resolution=ha.SIM_RESOLUTION)


if __name__ == '__main__':
    main()
