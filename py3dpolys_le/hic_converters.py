#! /usr/bin/env python

import argparse
import logging
import os
import sys

import cooler
import numpy as np
import pandas as pd

from . import hic_converters as hc
from ._version import __version__

# Initialization
logger = logging.getLogger(__name__)


def diff_matrix_to_cool(matrix1_file, matrix2_file, chr, resolution):
    bin_w = resolution
    if matrix1_file.endswith('gz'):
        matrix1_pd = pd.read_table(matrix1_file, delimiter='\t', index_col=0, comment='#', compression='gzip')
    else:
        matrix1_pd = pd.read_table(matrix1_file, delimiter='\t', index_col=0, comment='#')
    hic1 = matrix1_pd.to_numpy()
    np.fill_diagonal(hic1, 0)

    if matrix2_file.endswith('gz'):
        matrix2_pd = pd.read_table(matrix2_file, delimiter='\t', index_col=0, comment='#', compression='gzip')
    else:
        matrix2_pd = pd.read_table(matrix2_file, delimiter='\t', index_col=0, comment='#')
    hic2 = matrix2_pd.to_numpy()
    np.fill_diagonal(hic2, 0)

    hic_diff = hic2 - hic1
    cool_file = matrix1_file + f'_diff.{resolution}.cool'

    return hic_to_cool(hic_diff, chr, resolution, cool_file)


def matrix_to_cool(matrix_file, chr, resolution):
    if matrix_file.endswith('gz'):
        matrix_pd = pd.read_table(matrix_file, delimiter='\t', index_col=0, comment='#', compression='gzip')
    else:
        matrix_pd = pd.read_table(matrix_file, delimiter='\t', index_col=0, comment='#')

    logger.info(f"data.shape: ${matrix_pd.shape}")

    hic = matrix_pd.to_numpy()
    np.fill_diagonal(hic, 0)

    cool_file = matrix_file + f'.{resolution}.cool'
    return hic_to_cool(hic, chr, resolution, cool_file)


def hic_to_cool(hic, chr, resolution, cool_file):

    # build the cooler fields
    N = hic.shape[0]

    bins_index = [[chr, i * resolution, i * resolution + resolution] for i in range(N)]
    bins = pd.DataFrame(data=bins_index, columns=['chrom', 'start', 'end'])  # , dtype=np.dtype([('','','')]))

    pixels_bin1_id = []
    pixels_bin2_id = []
    pixels_count = []

    tot_iter = (N - 1) * N / 2
    iter = 0
    for bin1_id in range(N - 1):
        for bin2_id in range(bin1_id + 1, N):
            iter += 1
            progress = (iter / tot_iter) * 100
            if (progress % 10) == 0:
                logger.info(f'pixels progress: {progress}%')
            count = hic[bin1_id, bin2_id]
            if count != 0:
                # pixels_pd = pixels_pd.append({'bin1_id': np.int64(bin1_id), 'bin2_id': np.int64(bin2_id), 'count': count}, ignore_index=True)
                pixels_bin1_id.append(np.int64(bin1_id))
                pixels_bin2_id.append(np.int64(bin2_id))
                pixels_count.append(count)

    pixels_dic = {'bin1_id': pixels_bin1_id, 'bin2_id': pixels_bin2_id, 'count': pixels_count}
    metadata = {'format': 'HDF5::Cooler',
                'format-version': '0.8.6',
                'bin-type': 'fixed',
                'bin-size': resolution,
                'storage-mode': 'symmetric-upper',
                'genome-assembly': 'ce11',
                'generated-by': __name__ + '-' + __version__,
                # 'creation-date': datetime.date.today()
                }

    count_dtypes = {'count': 'float64'}
    cooler.create_cooler(cool_file, bins=bins, pixels=pixels_dic, dtypes=count_dtypes, ordered=True, metadata=metadata)
    return cool_file
    # problem with showing .cool file in higlass but with .mcool it works


def cool_to_matrix(cool_file, balanced=False):
    mat_cooler = cooler.Cooler(f'{cool_file}::/')
    chrom = mat_cooler.chromnames[-1]
    matrix_file = cool_file + f'_{chrom}.matrix.tsv'
    name = "wt_N2_Brejc2017_5000"
    species = "ce11"
    mat = mat_cooler.matrix(balance=balanced).fetch(mat_cooler.chromnames[0])  # TODO all chromosomes
    res = mat_cooler.binsize
    bin_names = [f'{name}|{species}|{chrom}:{bi*res}-{(bi+1)*res}' for bi in range(mat.shape[0])]
    mat_df = pd.DataFrame(data=mat, index=bin_names, columns=bin_names)
    result = mat_df.to_csv(matrix_file, sep='\t')

    return matrix_file


def matrix_to_mcool(matrix_file, chr, resolution, factors):
    cool_file = matrix_file + f'.{resolution}.cool'
    if not os.path.isfile(cool_file):
        matrix_to_cool(matrix_file, chr, resolution)  # == cool_file
    resolutions = [int(i * resolution) for i in factors]
    mcool_file = matrix_file + f'.{resolutions[0]}.mcool'
    cooler.zoomify_cooler(cool_file, mcool_file, resolutions=resolutions, chunksize=int(10e6))
    return mcool_file


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)

    p = argparse.ArgumentParser()
    p.add_argument("-cf", "--cool_file", default=f"{hc.PUBLISHED_FOLDER}/wt_N2_Brejc2017_5000.cool",
                   help="Experimental cool file to plot")
    p.add_argument("-o", "--output_folder", default=".", help="output folder")
    p.add_argument("-b", "--balanced", action='store_true', help="Balanced or not")
    p.add_argument("-z", "--z_score", action='store_true', help="Z-score normalized otherwise original Hi-C")
    p.add_argument("-chrs", "--chr_synonyms", nargs="+", default=['chrX', 'X', '6'],
                   help="List of chromosome synonyms to look for")
    p.add_argument("-c", "--cmap", default="hot_r ",
                   help="Color map: cool, hot_r, gist_heat_r, afmhot_r, YlOrRd, Greys, gist_yarg, seismic")
    p.add_argument("-f", "--file_extension", default="png", help="File format extension: png, tif")
    args = p.parse_args(sys.argv[1:])
