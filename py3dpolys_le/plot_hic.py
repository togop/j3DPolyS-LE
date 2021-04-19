#! /usr/bin/env python

import argparse
import logging
import os
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt
import glob

# Initialization
logger = logging.getLogger(__name__)


def run(output_folder, cmap, file_ext):
    logger.info(f'Plotting HiC for {output_folder} with color map: {cmap} in file format: {file_ext} ...')

    hic_files = glob.glob(os.path.join(output_folder, "hic*.hdf5"))
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

            plot_file = os.path.join(output_folder, os.path.basename(hic_file).replace('.hdf5', f'_{cmap}.{file_ext}'))

            fig.savefig(plot_file, dpi=200)
            plt.close()
            print(f'Hic plot saved in file {plot_file}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output_folder", default=".",
                   help="'Analysis' step output folder containing raw hic_*.hdf5 files.")
    p.add_argument("-c", "--cmap", default="hot_r",
                   help="Color map: cool, hot_r, gist_heat_r, afmhot_r, YlOrRd, Greys, gist_yarg")
    p.add_argument("-f", "--file_extension", default="png", help="Image Ffile format extension: png, tif")
    args = p.parse_args(sys.argv[1:])

    run(args.output_folder, args.cmap, args.file_extension)


if __name__ == '__main__':
    main()
