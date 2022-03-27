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
import dask.array as da
from scipy.signal import savgol_filter, find_peaks, find_peaks_cwt

# Initialization
logger = logging.getLogger(__name__)

DEFAULT_CMAP = "hot_r"
# tuned for simulation HiCs
DEFAULT_RESOLUTION = 2000
DEFAULT_HIC_WILDCARD = "hic*.hdf5"
DEFAULT_CLIM = [-2.75, 0]  # found to be mostly OK for our experimental and simulation HiCs
DEFAULT_PLOT_FORMAT = "png"
CHR_X_SYNONYMS = ['6', 'chrX', 'X']
DEFAULT_TITLE = "Hi-C for measurement {hic_file}"
DEFAULT_DPI = 150

# default simulation resolution
SIM_RESOLUTION = 2000

# experimental: extra optional
crop = None  # [9000//2, 15000//2]


def cli_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output_folder", default=".",
                   help="'Analysis' step output folder containing raw hic_*.hdf5 files.")
    p.add_argument("-w", "--hic_wildcard", default=DEFAULT_HIC_WILDCARD,
                   help="Hi-C files wildcard. Supported hic formats: *.hdf5 (simulation's format), *.cool, *.mcool")
    p.add_argument("--hic_chrs", nargs='*', default=CHR_X_SYNONYMS,
                   help="Synonyms of the chromosome from Hi-C matrixes to be ploted. "
                        "Default: synonyms for chrX in c.elegans: chrX, X, 6 ")
    p.add_argument("-r", "--resolution", default=DEFAULT_RESOLUTION, type=int, help="Hi-C data resolution in bp.")
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
            if resolution > SIM_RESOLUTION:  # need binning
                hic = hic_coarsen(hic, hic_file, resolution)
            res = SIM_RESOLUTION
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


def hic_coarsen(hic, hic_h5, resolution, from_resolution=SIM_RESOLUTION):
    factor = resolution // from_resolution
    if factor == 1:
        return hic   # no coarsen
    if resolution % from_resolution != 0:
        logger.error(f'Cannot do binning from resolution {from_resolution} to {resolution} '
                     f'for hic file {hic_h5}')
    # coarsen
    dhic = da.from_array(hic)  # , chunks=hic.shape[0] // factor)
    dhic_binned = da.coarsen(np.sum, dhic, {0: factor, 1: factor})
    hic = dhic_binned.compute()
    return hic


def run(output_folder, hic_wildcard=DEFAULT_HIC_WILDCARD, resolution=DEFAULT_RESOLUTION, balanced=False,
        hic_chrs=CHR_X_SYNONYMS, cmap=DEFAULT_CMAP, clim=DEFAULT_CLIM, title=DEFAULT_TITLE,
        plot_format=DEFAULT_PLOT_FORMAT):
    logger.info(f'Plotting HiC for {output_folder} with color map: {cmap} in file format: {plot_format} ...')
    if not cmap:
        cmap = DEFAULT_CMAP
    if not plot_format:
        plot_format = DEFAULT_PLOT_FORMAT

    hic_files = sorted(glob.glob(os.path.join(output_folder, hic_wildcard)))
    for hic_file in hic_files:
        hic = get_hic(hic_file, resolution, balanced, hic_chrs)
        print(f"data.shape: ${hic.shape}")

        fig = plt.figure()

        # crop desired section
        if crop:
            if len(crop) == 2:
                hic = hic[crop[0]:crop[1], crop[0]:crop[1]]
            else:
                hic = hic[crop[0]:, crop[0]:]

        hic_log = np.log10(hic)

        res_kb = resolution/1000
        plt.gca().get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: int(x*res_kb)))
        plt.gca().get_yaxis().set_major_formatter(FuncFormatter(lambda y, p: int(y*res_kb)))

        plt.imshow(hic_log, interpolation='nearest', cmap=cmap)  # color scale
        # cbar_h = plt.colorbar()

        # # START: loop anchors finder: for experimental HiCs
        # hic_aggregated = np.log10(hic.sum(axis=0)) * -1000
        # hic_aggregated[hic_aggregated == -np.inf] = 0
        # hic_aggregated[hic_aggregated == np.inf] = 0
        # hic_aggregated_smooth = savgol_filter(hic_aggregated, 7, 3) + 700
        # plt.plot(hic_aggregated_smooth)
        # # peaks, _ = find_peaks(1000-hic_aggregated_smooth, height=170, width=6)  # np.arange(10, 20))
        # # resolution: 20000
        # #peaks = find_peaks_cwt(-hic_aggregated_smooth, widths=np.arange(5, 30))
        # peaks = find_peaks_cwt(-hic_aggregated_smooth, widths=np.arange(5, 27))
        # plt.plot(peaks, hic_aggregated_smooth[peaks], "x")
        # # END: loop anchors finder

        print(plt.rcParams['axes.prop_cycle'].by_key()['color'])

        cbar_h = plt.colorbar()
        # mu = np.mean(hic)
        # sd = np.std(hic)
        # qt = np.quantile(hic, 0.50)
        # plt.clim(np.log10(mu - 1.64*sd), np.log10(mu + 1.64*sd))
        # plt.clim(np.log2(mu + 1.25*sd), 0)
        if clim and len(clim) > 1:
            plt.clim(clim[0], clim[1])
        hic_file_clean = hic_file.replace('_hic_003.hdf5', '').replace('./', '')   # extra clean
        plt.title(title.format(hic_file=hic_file_clean, output_folder=output_folder, resolution=resolution, balanced=balanced))

        # cbar_h.ax.tick_params(labelsize=11)
        # plt.show()

        basename = os.path.basename(hic_file)
        extension = basename.split(".")[-1]
        crop_ext = f'_c{crop[0]}{f"-{crop[1]}" if len(crop)==2 else ""}' if crop else ""
        plot_file = os.path.join(output_folder, basename.replace(f'.{extension}', f'{crop_ext}_{cmap}.{plot_format}'))

        fig.savefig(plot_file, dpi=DEFAULT_DPI)
        plt.close()
        print(f'Hic plot saved in file {plot_file}')


def main():
    args = cli_parser().parse_args(sys.argv[1:])

    run(args.output_folder, hic_wildcard=args.hic_wildcard, resolution=args.resolution, balanced=args.balanced,
        hic_chrs=args.hic_chrs, cmap=args.cmap, clim=args.clim, title=args.title, plot_format=args.plot_format)


if __name__ == '__main__':
    main()
