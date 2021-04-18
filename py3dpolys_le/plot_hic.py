#! /usr/bin/env python

import argparse
import logging
import os
import sys

import h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Initialization
logger = logging.getLogger(__name__)


def run(input_dat, output_folder, cmap, file_ext, m=None, chi2_alpha=1):
    input_dat_pd = pd.read_csv(input_dat, delim_whitespace=True, encoding='utf-8', names=['value', 'name'])
    #input_dat_pd.columns = ['value', 'name']
    input_dat_pd.set_index('name', inplace=True)

    Nmeas = int(input_dat_pd.loc['::Nmeas'][0])
    digits = int(np.log10(Nmeas)) + 1
    # I = 42

    logger.info(f'Plotting HiC for {Nmeas} measurements with color map: {cmap} in file format: {file_ext} ...')

    s_alpha = ''
    if m is not None:
        start_m = m
        end_m = m
    else:
        start_m = 1
        end_m = Nmeas

    sim_folder = '/'.join(os.path.abspath(output_folder).split('/')[-3:])
    for i in range(start_m, end_m + 1):
    #if i > 0:
        hic_file = os.path.join(output_folder, f'hic_{i:03}.hdf5')
        if not os.path.isfile(hic_file):  # legacy names
            hic_file = os.path.join(output_folder, f'hic_{i}.hdf5')
        with h5py.File(hic_file, 'r') as f:
            print(f"plotting HiC: {hic_file} , Keys: {f.keys()}")
            a_group_key = list(f.keys())[0]

            # Get the data
            data = list(f[a_group_key])

            hic = np.array(data)
            np.fill_diagonal(hic, 0)
            print(f"data.shape: ${hic.shape}")

            if chi2_alpha != 1:
                hic = hic * chi2_alpha
                s_alpha = f'_a{chi2_alpha:.4f}'

            fig = plt.figure()

            hic_log = np.log10(hic)

            plt.imshow(hic_log, interpolation='nearest', cmap=cmap)  # color scale
            cbar_h = plt.colorbar()

            print(plt.rcParams['axes.prop_cycle'].by_key()['color'])

            mu = np.mean(hic)
            sd = np.std(hic)
            #qt = np.quantile(hic, 0.50)
            #plt.clim(np.log10(mu - 1.64*sd), np.log10(mu + 1.64*sd))
            #plt.clim(np.log2(mu + 1.25*sd), 0)
            plt.clim(-2.75, 0)
            plt.title(f'Hi-C for measurement {i:0{digits}} from simulation:\n{sim_folder}')

            #cbar_h.ax.tick_params(labelsize=11)
            # plt.show()

            fig.savefig(os.path.join(output_folder, f'hic_{i:03}_{cmap}{s_alpha}.{file_ext}'), dpi=200)
            plt.close()

            # PCA
            # hic_z = StandardScaler().fit_transform(hic)
            #
            # pca = PCA(n_components=2)
            # pcs = pca.fit_transform(hic_z)
            # print(f"pcs.shape: ${pcs.shape}")
            #
            # plt.plot(range(pcs.shape[0]), pcs[:, 0])   #   1:hic.shape[1],
            # plt.show()
            #
            # print("pcs 2")
            # plt.plot(range(pcs.shape[0]), pcs[:, 1])
            #
            # print("end")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input_dat", default="../input.dat", help="Input.dat file used from a simulation.")
    p.add_argument("-o", "--output_folder", default=".",
                   help="'Analysis' step output folder containing raw hic_*.hdf5 files.")
    p.add_argument("-c", "--cmap", default="hot_r",
                   help="Color map: cool, hot_r, gist_heat_r, afmhot_r, YlOrRd, Greys, gist_yarg")
    p.add_argument("-f", "--file_extension", default="png", help="Image Ffile format extension: png, tif")
    p.add_argument("-m", "--measurement", default=None, help="Specific measurement to plot.", type=int)
    p.add_argument("-a", "--chi2_alpha", default=1, type=float,
                   help="In combination with the -m parameter to specify the chi2-min alpha coefficient "
                        "used as a multiplication factor, by which all Hi-C values will be multiplied.")
    args = p.parse_args(sys.argv[1:])

    run(args.input_dat, args.output_folder, args.cmap, args.file_extension,
        m=args.measurement, chi2_alpha=args.chi2_alpha)

    # /Volumes/imaging.data/ttgitchev/dcc-extrusion_sims/s72000_t100_Nlef600_brn5000m/out
    # /Volumes/ttgitchev/dcc-extrusion_sims/mex_s144000_t100_Nlef600_brn0/out

    # https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html

    #hic_file = 'hic_11.hdf5'
    #hic_coll_file = 'hic_41.cool'

    # hic2cool_convert(hic_file, hic_coll_file, 2000)
    # h5file = h5py.File(hic_coll_file, 'r')
    # ### will give you the cooler object with resolution = 10000 bp
    # hic_coll = cooler.Cooler(h5file)


#  hicPlotMatrix -m hic_41.hdf5.2000.cool  -o hic_41.hdf5.2000.png --log1p     # --region 1:20000000-80000000
#  hicFindTADs -m hic_41.hdf5.2000.cool --outPrefix hic_41 --correctForMultipleTesting None --minBoundaryDistance 4000 --numberOfProcessors 1
#  hicPCA -m hic_41.hdf5.2000.cool -o hic_41_pca1.bw hic_41_pca2.bw --format bigwig
#  hicPlotMatrix -m hic_41.hdf5.2000.cool -o hic_41_pca1.png --perChr --bigwig hic_41_pca1.bw

if __name__ == '__main__':
    main()
