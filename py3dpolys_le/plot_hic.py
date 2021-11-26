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
import cooler

# Initialization
logger = logging.getLogger(__name__)

DEFAULT_CMAP = "hot_r"
DEFAULT_CLIM = [-2.75, 0]  # tuned for simulation HiC
DEFAULT_PLOT_FORMAT = "png"
CHR_X_SYNONYMS = ['6', 'chrX', 'X']
DEFAULT_TITLE = "Hi-C for measurement {hic_file}"  # as used for simulation HiC


def cli_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output_folder", default=".",
                   help="'Analysis' step output folder containing raw hic_*.hdf5 files.")
    p.add_argument("-w", "--files_wildcard", default="hic*.hdf5", help="Hi-C files wildcard.")
    p.add_argument("--hic_chrs", nargs='*', default=CHR_X_SYNONYMS,
                   help="Synonyms of the chromosome from Hi-C matrixes to be ploted. "
                        "Default: synonyms for chrX in c.elegans: chrX, X, 6 ")
    p.add_argument("-r", "--resolution", default=2000, type=int, help="Hi-C data resolution in bp.")
    p.add_argument("-b", "--balanced", action='store_true', help="Get balanced if available, for .cool and .mcool.")
    p.add_argument("-c", "--cmap", default=DEFAULT_CMAP,
                   help="Color map: cool, hot_r, gist_heat_r, afmhot_r, YlOrRd, Greys, gist_yarg. "
                        f"Default: {DEFAULT_CMAP}")
    p.add_argument("--clim", nargs='*', type=float, default=DEFAULT_CLIM,
                   help="Set the color limits of the current image, see matplotlib.pyplot.clim parameter."
                        f"Default: {DEFAULT_CLIM[0]}, {DEFAULT_CLIM[1]}")
    p.add_argument("-f", "--plot_format", default=DEFAULT_PLOT_FORMAT,
                   help=f"Image file format extension: png, tif, svg. Default: {DEFAULT_PLOT_FORMAT}")
    p.add_argument("--title", default=DEFAULT_TITLE,
                   help="Plot title. It could be template containing parameters in {}: hic_file, output_folder, "
                        "resolution, balanced"
                        f". Default: {DEFAULT_TITLE}")
    return p


def get_hic(hic_file, resolution, balanced, hic_chrs):
    if hic_file.endswith('.hdf5'):
        # directly from .hdf5
        with h5py.File(hic_file, 'r') as f:
            a_group_key = list(f.keys())[0]
            logger.info(f"get_hic for {hic_file} HDF5 Keys {f.keys()} and dimentions {a_group_key}")
            # Get the data
            data = list(f[a_group_key])
            hic = np.array(data)
            np.fill_diagonal(hic, 0)
    elif hic_file.endswith('.cool') or hic_file.endswith('.mcool'):
        # use the root cool
        cooler_ref = f'{hic_file}::/' if hic_file.endswith('.cool') else f'{hic_file}::/resolutions/{resolution}'
        hic_cooler = cooler.Cooler(cooler_ref)
        hic_chr = list(set(hic_cooler.chromnames) & set(hic_chrs))[0]
        balance = (hic_cooler.bins()['weights'] is not None) if balanced else balanced
        hic = hic_cooler.matrix(balance=balance).fetch(hic_chr)
        if balanced:
            hic = np.nan_to_num(hic)  # nan -> 0
        np.fill_diagonal(hic, 0)
    return hic


def run(output_folder, files_wildcard, resolution, balanced, hic_chrs, cmap, clim, title, plot_format):
    logger.info(f'Plotting HiC for {output_folder} with color map: {cmap} in file format: {plot_format} ...')
    if not cmap:
        cmap = DEFAULT_CMAP
    if not plot_format:
        plot_format = DEFAULT_PLOT_FORMAT

    hic_files = sorted(glob.glob(os.path.join(output_folder, files_wildcard)))
    for hic_file in hic_files:
        hic = get_hic(hic_file, resolution, balanced, hic_chrs)
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
        if clim and len(clim) > 1:
            plt.clim(clim[0], clim[1])
        plt.title(title.format(hic_file=hic_file, output_folder=output_folder, resolution=resolution, balanced=balanced))

        # cbar_h.ax.tick_params(labelsize=11)
        # plt.show()

        basename = os.path.basename(hic_file)
        extension = basename.split(".")[-1]
        plot_file = os.path.join(output_folder, basename.replace(f'.{extension}', f'_{cmap}.{plot_format}'))

        fig.savefig(plot_file, dpi=200)
        plt.close()
        print(f'Hic plot saved in file {plot_file}')


def main():
    args = cli_parser().parse_args(sys.argv[1:])

    run(args.output_folder, args.files_wildcard, args.resolution, args.balanced, args.hic_chrs, args.cmap, args.clim,
        args.title, args.plot_format)


if __name__ == '__main__':
    main()
