#! /usr/bin/env python

import argparse
import logging
import sys
import os
import numpy as np
#import re
#from scipy.stats import spearmanr, pearsonr

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

PLOT_3D = '3d'
PLOT_HMAP = 'hmap'

# Initialization
SUBSTR_UNIDIR = '_unidir'
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("").setLevel(logging.INFO)
logger = logging.getLogger(__name__)

p = argparse.ArgumentParser()
p.add_argument("-f", "--stats_file", default="./sim_stats.csv",
               help="Simulation statistics’ repository file (sim_stats.csv).")
p.add_argument("-o", "--output_folder", default=".",
               help="Output folder to save plots.")
p.add_argument("-e", "--file_extension", default="png", help="Image file format extension: png, tif, svg.")
p.add_argument("-z", "--z_column", default="chi2_log",
               help="Column from a --stats_file to be plotted on the z-axis.")
# p.add_argument("-u", "--unidirectional", help="Show only unidirectional mode for LEFs.",
#               action='store_true')
# p.add_argument("-b", "--bidirectional", help="Show only bidirectional mode for LEFs.",
#               action='store_true')
p.add_argument("-p", "--plot_mode", choices=[PLOT_3D, PLOT_HMAP], default=PLOT_HMAP,
               help="Plot method. 3d - 3-dimensional scatter plot; hmap - 2D heatmap plot.")
p.add_argument("-c", "--cmap", default="YlGnBu_r",   # "YlGnBu_r" appropriate for 2D, "cool" for heatmap
               help="Color map: YlGnBu, cool, hot_r, gist_heat_r, afmhot_r, YlOrRd, Greys, gist_yarg")
# p.add_argument("--list_km", nargs='+', default=, help="Shows data only for list of LEF velocities (km)")
p.add_argument("--list_nlef", nargs='*', default=[], help="Show data only for the list of of LEFs occupancy (Nlefs).")
p.add_argument("-lr", "--list_contact_radii", nargs='+', default=['2.84', '3.55'],
               help="Show data only for the list of contact radii with a different color for every contact radius.")
args = p.parse_args(sys.argv[1:])

logger.info(f'start with parameters: {args}')

# see https://matplotlib.org/3.1.3/tutorials/colors/colormaps.html  for palette information
CONTACT_RADIUS_PALETTE = plt.get_cmap('tab10')


def plot_for_radius(subplot, r, z_col, color, list_nlef=None):
    # global sim_stats_pd, x, y, z_chi2
    # TODO find better way instead of using fixed strings: z-loop, unidir
    sim_stats_pd = sim_stats_all_pd[(sim_stats_all_pd.radius_contact == r)]

    if list_nlef:
        list_nlef_i = [int(n) for n in list_nlef]
        sim_stats_pd = sim_stats_pd[sim_stats_pd.nlef.isin(list_nlef_i)]

    sim_stats_pd = sim_stats_pd.sort_values(['sim_hic_file'])

    # sim_stats_pd_2 = sim_stats_all_pd[sim_stats_all_pd_2.radius_contact == r].sort_values(['sim_hic_file'])
    x = sim_stats_pd[['nlef']]
    y = sim_stats_pd[['km']]
    z = sim_stats_pd[[z_col]]
    # z_val = np.array(z.values)   # - np.array(sim_stats_pd_2[z_lab].values)
    if args.plot_mode == PLOT_HMAP:
        heatmap_data = pd.pivot_table(sim_stats_pd, values=z_col,
                                       index=['nlef'],
                                       columns='km')

        heatmap_data_tr = heatmap_data.rename(columns=lambda km: f'{km:.2e}')
        # heatmap_data_tr = heatmap_data.transform(axis=1, func=lambda km: f'{km:.2f}')

        return sns.heatmap(heatmap_data_tr, annot=True, fmt='.2f', cmap=args.cmap, cbar_kws={'label': z_col})
    else:
        return subplot.scatter(x, y, z, color=color, marker='o')
        # subplot.plot_surface(x, y, z, color=color) #, rstride=1, cstride=1, cmap='viridis', edgecolor='none')


def plot_3d_stats(z_col, list_nlef, subplot):
    # fig = plt.figure()
    # plt_1 = fig.add_subplot(111, projection='3d')
    legend = []
    color_shift = 0
    len_radii = len(args.list_contact_radii)
    loop_radii(args.list_contact_radii, list_nlef, color_shift, legend, subplot, z_col)

    if args.plot_mode == PLOT_3D:
        subplot.set_xlabel('Nlef')
        subplot.set_ylabel('km')
        subplot.set_zlabel(z_col)
        subplot.legend(legend)
    # plt.show()


def loop_radii(radii, list_nlef, color_shift, legend, subplot, z_col):
    for i, r in enumerate(radii):
        print(f'plot radius:{r} with color:{color_shift + i}')
        if args.plot_mode == PLOT_HMAP:
            fig = plt.figure()
        subplot = plot_for_radius(subplot, r=r, z_col=z_col, color=CONTACT_RADIUS_PALETTE(color_shift + i)
                                  , list_nlef=list_nlef)

        if args.plot_mode == PLOT_HMAP:
            basename = os.path.basename(args.stats_file)
            fig.suptitle(f'Heatmap of \'{z_col}\' with contact radius:{r} for: \n{basename}')
            subplot.set_xlabel('km')
            subplot.set_ylabel('Nlef')

            # plt.show()

            file_name: str = f'{basename}_{args.plot_mode}_{z_col}_r{r}.{args.file_extension}'
            fig.savefig(os.path.join(args.output_folder, file_name), dpi=200)
            logger.info(f'Plot saved in {file_name}')
            plt.close()


sim_stats_all_pd = pd.read_csv(args.stats_file, delimiter=',', encoding='utf-8', header=0,
                               dtype={'radius_contact': str})
# if you want to compare:
# sim_stats_all_pd_2 = pd.read_csv('./sim_stats_lyon_chi2-sem.csv', delimiter=',', encoding='utf-8', header=0)

if args.plot_mode == PLOT_3D:
    fig = plt.figure()
    # plt_1 = fig.add_subplot(221, projection='3d')

    plt_1 = fig.add_subplot(111, projection='3d')
else:
    plt_1 = None   # later one for every plot
#    plt_1 = fig.add_subplot(111)

# z_cols: chi2_log, alpha_log, chi2_log_1tad, alpha_log_1tad, chi2_lin_1tad, alpha_lin_1tad
plot_3d_stats(args.z_column, args.list_nlef, plt_1)

# plt_2 = fig.add_subplot(122, projection='3d')
# plot_3d_stats('chi2_log', plt_2)

# TODO make build help CLI; clean code


def main():
    # plt.show()  # only for interactive mode to choose the right perspective

    file_name = f'{os.path.basename(args.stats_file)}_{args.plot_mode}_{args.z_column}' \
                f'_r{"_r".join(args.list_contact_radii)}.{args.file_extension}'
#                f'{"_u" if args.unidirectional else ""}{"_b" if args.bidirectional else ""}' \
    fig.savefig(os.path.join(args.output_folder, file_name), dpi=200)

    logger.info(f'Plot saved in {file_name}')


if args.plot_mode == PLOT_3D:
    main()
