#! /usr/bin/env python

import argparse
import logging
import os
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
import glob

# Initialization
logger = logging.getLogger(__name__)

DEFAULT_CMAP = "hot_r"
DEFAULT_PLOT_FORMAT = "png"


def cli_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output_folder", default=".",
                   help="'Analysis' step output folder containing raw hic_*.hdf5 files.")
    p.add_argument("-r", "--resolution", default=2000, help="Hi-C data resolution in bp.")
    p.add_argument("-c", "--cmap", default=DEFAULT_CMAP,
                   help="Color map: cool, hot_r, gist_heat_r, afmhot_r, YlOrRd, Greys, gist_yarg. "
                        f"Default: {DEFAULT_CMAP}")
    p.add_argument("-f", "--plot_format", default=DEFAULT_PLOT_FORMAT,
                   help=f"Image file format extension: png, tif, svg. Default: {DEFAULT_PLOT_FORMAT}")
    return p


def run(output_folder, resolution, cmap, plot_format):
    logger.info(f'Plotting HiC for {output_folder} with color map: {cmap} in file format: {plot_format} ...')
    if not cmap:
        cmap = DEFAULT_CMAP
    if not plot_format:
        plot_format = DEFAULT_PLOT_FORMAT

    hic_files = sorted(glob.glob(os.path.join(output_folder, "hic*.hdf5")))
    for hic_file in hic_files:
        with h5py.File(hic_file, 'r') as f:
            print(f"plotting HiC: {hic_file} , Keys: {f.keys()}")
            a_group_key = list(f.keys())[0]

            # Get the data
            data = list(f[a_group_key])

            hic = np.array(data)
            np.fill_diagonal(hic, 0)
            print(f"data.shape: ${hic.shape}")

            fig = plt.figure()

            hic_log = np.log10(hic)

            res_kb = resolution//1000
            plt.gca().get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: int(x*res_kb)))
            plt.gca().get_yaxis().set_major_formatter(FuncFormatter(lambda y, p: int(y*res_kb)))

            plt.imshow(hic_log, interpolation='nearest', cmap=cmap)  # color scale
            # cbar_h = plt.colorbar()

            print(plt.rcParams['axes.prop_cycle'].by_key()['color'])

            # mu = np.mean(hic)
            # sd = np.std(hic)
            # qt = np.quantile(hic, 0.50)
            # plt.clim(np.log10(mu - 1.64*sd), np.log10(mu + 1.64*sd))
            # plt.clim(np.log2(mu + 1.25*sd), 0)
            plt.clim(-2.75, 0)
            plt.title(f'Hi-C for measurement {hic_file} from simulation:\n{output_folder}')

            # cbar_h.ax.tick_params(labelsize=11)
            # plt.show()

            plot_file = os.path.join(output_folder, os.path.basename(hic_file).replace('.hdf5', f'_{cmap}.{plot_format}'))

            fig.savefig(plot_file, dpi=200)
            plt.close()
            print(f'Hic plot saved in file {plot_file}')


def main():
    args = cli_parser().parse_args(sys.argv[1:])

    run(args.output_folder, args.resolution, args.cmap, args.plot_format)


if __name__ == '__main__':
    main()
