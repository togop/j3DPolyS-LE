#! /usr/bin/env python

import argparse
import logging
import os
import re

import cooler
import h5py
import numpy as np
import pandas as pd
import scipy.io
import sys

from py3dpolys_le import _version

# Initialization
logger = logging.getLogger(_version.__name__)


def cli_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input_file",
                   help="input file name. Supported file types: .hdf5 (3DPolyS_LE format), .mat (MATLAB 2D matrix)")
    p.add_argument("-o", "--output_file", help="output file name. Allowed file types: .cool, .mcool")
    p.add_argument("-chr", help="Chromosome name to be used for. Default: chrS", default="chrS")
    p.add_argument("-r", "--resolutions", nargs='+', default=[2000], type=int,
                   help="List of resolutions for .mcool output file. "
                        "Default: 2000, for a .cool file in a simulation's resolution.")
    # p.add_argument("-f", "--foctors", help="Chromosome name to be used for")
    return p


def read_hic_hdf5(hic_hdf5):
    with h5py.File(hic_hdf5, 'r') as f:
        # List all groups
        logger.info(f"Read hic file {hic_hdf5} with Keys: {f.keys()} ...")
        a_group_key = list(f.keys())[0]
        data = list(f[a_group_key])
        hic = np.array(data)
        return hic


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
                'generated-by': _version.__name__ + '-' + _version.__version__,
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


def hdf5_to_cooler():
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)

    p = cli_parser()
    args = p.parse_args(sys.argv[1:])

    hic = read_hic_hdf5(args.input_file)
    if args.output_file.endswith('.cool'):
        hic_to_cool(hic=hic, chr=args.chr, resolution=args.resolutions[0], cool_file=args.output_file)
    elif args.output_file.endswith('.mcool'):
        hic_to_mcool(hic, chr=args.chr, resolutions=args.resolutions, mcool_file=args.output_file)


def hic_to_mcool(hic, chr, resolutions, mcool_file):
    cool_file = re.sub(r'\.mcool', '.cool', mcool_file)  # hic to compare with
    hic_to_cool(hic, chr=chr, resolution=resolutions[0], cool_file=cool_file)
    # using CLI
    res_str = str(resolutions).strip('[]')
    cmd = f'cooler zoomify -o {mcool_file} -c 10000000 -r \'{res_str}\' {cool_file}'
    logging.info(f'call: {cmd}')
    os.system(cmd)
    # using API: there was some problems sometime
    # try:
    #    cooler.zoomify_cooler(cool_file, args.output_file, resolutions=args.resolutions, chunksize=int(10e6))
    # except:
    #    logger.warning(f'Failed to zoomify file {cool_file}! Will try CLI cooler zoomify ...')


def mat_to_cooler():
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)

    p = cli_parser()
    args = p.parse_args(sys.argv[1:])

    hic = scipy.io.loadmat(args.input_file).get('m')
    # hic_to_cool(hic=hic, chr='chr1', resolution=2000, cool_file='/Users/todor/unibe/git/3DPolyS-LE/test/data/fakematrix_TAD300_600k_peaked.mcool')
    if args.output_file.endswith('.cool'):
        hic_to_cool(hic=hic, chr=args.chr, resolution=args.resolutions[0], cool_file=args.output_file)
    elif args.output_file.endswith('.mcool'):
        hic_to_mcool(hic, chr=args.chr, resolutions=args.resolutions, mcool_file=args.output_file)
    print("DONE")


def main():
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)

    p = cli_parser()
    args = p.parse_args(sys.argv[1:])
    if args.input_file.endswith('.hdf5'):
        hic = read_hic_hdf5(args.input_file)
    elif args.input_file.endswith('.mat'):
        hic = scipy.io.loadmat(args.input_file).get('m')
    else:
        logging.error(f'Unsupported input format {args.input_file}. '
                      f'Supported file types: .hdf5 (3DPolyS_LE format), .mat (MATLAB 2D matrix)')
        exit(1)

    if args.output_file.endswith('.cool'):
        hic_to_cool(hic=hic, chr=args.chr, resolution=args.resolutions[0], cool_file=args.output_file)
    elif args.output_file.endswith('.mcool'):
        hic_to_mcool(hic, chr=args.chr, resolutions=args.resolutions, mcool_file=args.output_file)
    elif args.output_file.endswith('.mat'):
        mdic = {"m": hic}
        scipy.io.savemat(args.output_file, mdic)
    else:
        logging.error(f'Unsupported output format {args.output_file}. '
                      f'Supported file types: .cool, .mcool, .mat')
        exit(1)


if __name__ == "__main__":
    main()
    # hdf5_to_cooler()
