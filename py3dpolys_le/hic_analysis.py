import bisect
import copy
import csv
import fnmatch
import logging
import math
import os
import pathlib
import re
import subprocess

import cooler
import dask.array as da
import h5py
import numpy as np
import pandas as pd
import pyranges as pr
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy import stats

# from matplotlib import cm
# from matplotlib.colors import LinearSegmentedColormap

# Initialization
logger = logging.getLogger(__name__)

SIM_RESOLUTION = 2000
EXP_RESOLUTION = 10000

RESOLUTION = 10000

RES_FACTORS = [1, 2, 4]  # , 8, 16]

# CMAP = 'hot'
CMAP = 'YlOrRd'

# https://matplotlib.org/3.1.3/tutorials/colors/colormaps.html
DECAY_CMAP = 'tab20'  # ['spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia', 'pink', 'cooper']
DECAY_PALETTE = plt.get_cmap(DECAY_CMAP)
DECAY_PALETTE_SHIFT = 2
DECAY_COLOR_MAX = 20
DECAY_PLOT_ALPHA = 0.5
DECAY_USE_HDF5 = True
DECAY_USE_HDF5_COOL = False   # use the root cool file(produced by the hdf5): used only for demonstration purpose

SAVE_DECAY_PROBABILITY = True

CHI2_MODE_LOG = 'log'  # logarithmic so that #10^x == 100
CHI2_MODE_LINEAR = 'linear'  # linear: 20*[1..100]
CHI2_MODE_LOG_ZOOM = False

CHI2_RANGE_START = 20   # in SIM_RESOLUTION
CHI2_RANGE_END = 2000
CHI2_RANGE_NUM = 100

CHI2_USE_SEM = True
CHI2_USE_BALANCED = False  # True used only for plotting mixed comparisons

# PLOT_FORMAT = 'svg'
PLOT_FORMAT = 'png'

# simulation output files
DR_OUT = 'dr.out'
CONTACT_OUT = 'contact.out'
CONFIG_OUT = 'config.out'
NLEF_OUT = 'Nlef.out'
PROCESS_OUT = 'process.out'
# analyse output file
CHIP_OUT = 'chip_lef.out'
CHIP_BED_GRAPH = 'chip_lef.bedGraph'
XYZCONFIG_OUT = 'xyzconfig_001.out'   # at least the first one

ALPHA = r'$\alpha$'

# MC_HiC pipeline way
# clr_map = [cm.hot(ci) for ci in np.linspace(1.0, 0.0, 10)]
# CMAP = LinearSegmentedColormap.from_list('test', clr_map, N=10)

SIM_CHR_SYNONYMS = ['6', 'chrX', 'X', '1', 'chrI', 'I']
SIM_CHR = 'chrX'

CHR_X_SYNONYMS = ['6', 'chrX', 'X']
chr_x_size = 17718942
CHR_A_SYNONYMS = ['1', 'chrI', 'I']
chr_a_size = 15072434

CHR_SYNONYMS = CHR_X_SYNONYMS
chr_size = chr_x_size

DEFAULT_CHIP_CORRELATION = 'spearmanr'

PLOTS_FOLDER = 'plots'
DEMO_CLIM = [-1, 3]   # -2.75, 0]  #
DEFAULT_CLIM = DEMO_CLIM  # just from the DEMO: TODO set DEFAULT_CLIM=None
PLOT_COMP_TADS = False


def get_last_hic(output_folder: str):
    hic_h5s = fnmatch.filter(os.listdir(output_folder), 'hic_*.hdf5')
    if len(hic_h5s) == 0:
        logger.error(f'Not found hic_*.hdf5 files in {output_folder} found: {os.listdir(output_folder)}')
        exit(1)
    sim_hic_snapshot = np.sort(hic_h5s)[-1]  # 'hic_81.hdf5'  # test
    return os.path.join(output_folder, sim_hic_snapshot)


def hic_to_mcool(hic_file, chr, resolution, factors):
    cool_file = hic_file + f'.{resolution}.cool'
    if not os.path.isfile(cool_file):
        cool_file = hic_to_cooler(hic_file, chr, resolution)  # == cool_file
    resolutions = [int(i * resolution) for i in factors]
    mcool_file = hic_file + f'.{resolutions[0]}.mcool'
    try:
        cooler.zoomify_cooler(cool_file, mcool_file, resolutions=resolutions, chunksize=int(10e6))
    except:
        logger.warning(f'Problem to zoomify file {cool_file} so will regenerate it!')
        os.remove(cool_file)
        cool_file = hic_to_cooler(hic_file, chr, resolution)
        cooler.zoomify_cooler(cool_file, mcool_file, resolutions=resolutions, chunksize=int(10e6))
    return mcool_file


def balance_mcool(cool_file, resolutions, mcool_file):
    # cooler.zoomify_cooler(cool_file, mcool_file, resolutions=resolutions, chunksize=int(10e6))

    res_str = str(resolutions).strip('[]')
    cmd = f'cooler zoomify --balance --balance-args \'--convergence-policy store_nan\' -o {mcool_file} -c 10000000 -r \'{res_str}\' {cool_file}'
    os.system(cmd)

    logger.info(f' check balance for {mcool_file}')
    for res in resolutions:
        hic_cooler = cooler.Cooler(f'{mcool_file}::/resolutions/{res}')
        # cooler.balance_cooler(hic_cooler, map=mp.Pool(6).map)  # , ignore_diags=3, min_nnz=10)
        chr = list(set(hic_cooler.chromnames) & set(CHR_SYNONYMS))[0]
        hic_mat = hic_cooler.matrix(balance=True).fetch(chr)
        logger.info(f'resolution: {res} hic_mat.shape: {hic_mat.shape}')
    return mcool_file


def print_hdf5_structure(hdf5_file):
    logger.info(f'HDF5 structure of file ${hdf5_file}:')
    with h5py.File(hdf5_file, 'r') as f:
        # List all groups
        logger.info(f'Keys: ${f.keys()}')
        group_keys = list(f.keys())
        for key in list(group_keys):
            group = f[key]
            logger.info(f'group ${group}')
            if hasattr(group, 'keys'):
                logger.info(f'\tkeys: ${group.keys()}')
                for skey in list(group.keys()):
                    sgroup = group[skey]
                    logger.info(f'\tsub group ${sgroup}')
                    if hasattr(sgroup, 'keys'):
                        logger.info(f'\t\tkeys: ${sgroup.keys()}')


# print_hdf5_structure('/Users/todor/UniBern/master_project/MC-HiC-guppy/target/frag_files/frg_20190501_HIC6_7_barcode08_pass_WS235.hdf5')

def remove_duplicates(list, max=None):
    final_list = []
    for el in list:
        if (el <= max if max else True) and el not in final_list:
            final_list.append(el)
    return final_list


def hic_to_cooler(hic_file, chr=SIM_CHR, resolution=SIM_RESOLUTION):
    cool_file = f'{hic_file}.{resolution}.cool'
    with h5py.File(hic_file, 'r') as f:
        # List all groups
        logger.info(f"Converting hic file {hic_file} to {cool_file} with Keys: {f.keys()} ...")
        a_group_key = list(f.keys())[0]

        logger.info("a_group_key0: %s" % a_group_key)

        # Get the data
        data = list(f[a_group_key])

        hic = np.array(data)

        # np.fill_diagonal(hic, 0)  # no need: should be done already

        logger.info(f"data.shape: ${hic.shape}")

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
                    'assembly': 'ce11',
                    'generated-by': __name__,
                    # 'creation-date': datetime.date.today()
                    }

        count_dtypes = {'count': 'float64'}
        if not os.path.isfile(cool_file):
            cooler.create_cooler(cool_file, bins=bins, pixels=pixels_dic, dtypes=count_dtypes, ordered=True,
                                 metadata=metadata)
        else:
            logger.warning(f'Cooler file {cool_file} already exist and will not be replaced!')
        return cool_file
        # problem with showing .cool file in higlass but with .mcool it works


def save_chromosome_sizes(chr_lst, chr_size, file_name):
    chroms = pd.DataFrame(columns=['name', 'length'])
    with open(file_name, "w") as chromsizes_file:
        writer = csv.writer(chromsizes_file, delimiter='\t')
        # record_file.write("name    length\n")
        for chr_ind, chr_name in enumerate(chr_lst):
            chroms.loc[chr_ind] = [chr_name, chr_size[chr_ind]]
            writer.writerow([str(chr_ind + 1), chr_size[chr_ind]])
            # chromsizes_file.write(str(chr_name)+"   "+str(chr_size[chr_ind])+"\n")


def average_contact_prob(prob_mat, dist, plot=False, ax=plt):
    count = 0
    probs = []
    x = []
    y = []
    for i in range(dist, prob_mat.shape[0]):
        j = i - dist
        count += 1
        probs.append(prob_mat[i, j])
        if plot:
            x.append(i)
            y.append(j)

    # probs_filtered = list(filter(lambda x: x > 0, probs)) if count > 0 and filtered else probs
    # count = len(probs_filtered)
    avrg_prob = np.mean(probs, dtype=np.longdouble) if count > 0 else 0
    if CHI2_USE_SEM:
        sd_sem = stats.sem(probs) if count > 1 else 0
    else:
        sd_sem = np.std(probs, dtype=np.longdouble) if count > 1 else 0  # np.std(probs)**2 == np.var(probs)
    if plot:
        ax.plot(x, y)
    # avrg_prob_log = -np.log(max(avrg_prob, 1.e-8)) if avrg_prob < 1. else 0
    # std_prob_log = 0.1 * avrg_prob_log
    # return avrg_prob_log, std_prob_log  # , sem
    # float128: need highest precision
    return np.longdouble(avrg_prob), np.longdouble(sd_sem)  # std_prob  # , sem


def get_chi2_dist_range(chi2_mode, res, max=None):
    if chi2_mode == CHI2_MODE_LINEAR:
        linear_end = min(CHI2_RANGE_END, int(round(max * res / SIM_RESOLUTION))) if max else CHI2_RANGE_END
        dist_range = np.linspace(CHI2_RANGE_START, linear_end, num=CHI2_RANGE_NUM)  #
    else:  # chi2_mode == CHI2_MODE_LOG
        dist_range = CHI2_RANGE_START * np.logspace(0, np.log10(CHI2_RANGE_END/CHI2_RANGE_START),
                                                    num=CHI2_RANGE_NUM)  # // 10)  # if wanted less points
        # OLD WAY [(log_base ** x) for x in range(1, max_x + 1)]  # math.floor(math.exp(x))  # log_base**x
    dist_range = [int(round(x)) for x in (dist_range * SIM_RESOLUTION / res)]
    dist_range = np.asarray(remove_duplicates(dist_range, max))  # - 1  # TODO check also: filter by max; index start from 0(-1)
    return dist_range


def chi2_minimization(hic1_mat, cmp_hic_cooler, chrs, res, tads, comp_filename, plots_folder, norm=True,
                      chi2_mode=CHI2_MODE_LINEAR, hic_balance=CHI2_USE_BALANCED,
                      plot_cmap=CMAP, plot_format=PLOT_FORMAT):
    h1_title = 'simulation HiC'  # hic1cooler.info
    h2_title = 'experimental HiC'  # hic1cooler.info

    comp_chr = list(set(cmp_hic_cooler.chromnames) & set(chrs))[0]

    # total sums in the chi2_min
    tot_PS = 0
    tot_FS = 0
    tot_PFS = 0
    tads_chi2_min_df = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Value'])

    norm_term = 0
    for tadi in tads:
        tad_start = tadi[0]
        tad_end = tadi[1]
        tadi_size = (tad_end - tad_start) // res
        toti_PFS = 0
        toti_PS = 0
        toti_FS = 0
        if tadi_size == 0:
            logger.info(f'skip TAD:{tad_start}-{tad_end}, size:{tadi_size}')
        else:
            logger.info(f'calculate TAD:{tad_start}-{tad_end}, size:{tadi_size}')
            balanced = (cmp_hic_cooler.bins()['weights'] is not None) if hic_balance else hic_balance
            tadi_mat2 = (cmp_hic_cooler.matrix(balance=balanced).fetch((comp_chr, tad_start, tad_end - res)))  # -1 bead: fetch do inclusive start, end
            if balanced:
                tadi_mat2 = np.nan_to_num(tadi_mat2)
            tadi_mat1_start = tad_start // res
            tadi_mat1_end = tad_end // res
            tadi_mat1 = hic1_mat[tadi_mat1_start:tadi_mat1_end, tadi_mat1_start:tadi_mat1_end]

            if plots_folder:
                fig, (ax1, ax2) = plt.subplots(1, 2)
                try:
                    ax1.imshow(np.log10(tadi_mat1), cmap=plot_cmap, interpolation='nearest')
                    ax2.imshow(np.log10(tadi_mat2), cmap=plot_cmap, interpolation='nearest')
                except:
                    logger.warning(f'Problem to plot TAD: {tad_start}-{tad_end}')
            else:
                (ax1, ax2) = (0, 0)

            # log_base = 2  # 10 is too sparse  # e = math.exp(1)
            # max_x = max(1, np.int64(math.log(tadi_size, log_base)))
            # np.log((tad_end - tad_start) / res)  # /np.log(log_base))
            dist_range = get_chi2_dist_range(chi2_mode, res, tadi_size)
            # logger.info('chi2 used distances: ' + ', '.join(map(str, dist_range)))
            # dist_range = remove_duplicates(dist_range)   # so far not needed

            # norm_term += len(dist_range)
            if len(dist_range) < 1:
                logger.warning(
                    f'Empty chi2-dist-range for tad with size chi2(mode:{chi2_mode}, resolution:{res}), TADi_size:{tadi_size}[{tad_start}-{tad_end}]')
            tadi_norm_term = 0
            for dist in dist_range:

                (p_i, p_i_sdsem) = average_contact_prob(tadi_mat1, dist, plots_folder, ax1)
                (f_i, f_i_p_i_sdsem) = average_contact_prob(tadi_mat2, dist)
                if plots_folder:
                    average_contact_prob(tadi_mat2, dist, plots_folder, ax2)  # use just to plot
                sigma_i_2 = f_i_p_i_sdsem ** 2  # can use SD or SEM!!
                if sigma_i_2 != 0:
                    toti_PFS += (p_i * f_i) / sigma_i_2
                    toti_PS += (p_i ** 2) / sigma_i_2
                    toti_FS += (f_i ** 2) / sigma_i_2
                    norm_term += 1
                    tadi_norm_term += 1
                    logger.info(
                        f'chi2(dis:{dist}) : p_i: {p_i} * f_i: {f_i} = {p_i * f_i}, sigma_i^2= {sigma_i_2}, '
                        f'tot_FS={toti_FS}, tot_PFS={toti_PFS}, tot_PS={toti_PS}')
                else:
                    logger.info(
                        f'Skip chi2(dis:{dist}=tad_mats.shape:{tadi_mat2.shape[0]}) NAN sigma_i^2= {sigma_i_2},')
            tot_PFS += toti_PFS  # / max_x if norm else toti_PFS   # different way of normalization need to discuss
            tot_PS += toti_PS  # / max_x if norm else toti_PS
            tot_FS += toti_FS  # / max_x if norm else toti_FS

        tadi_chi2_min = (toti_FS - (toti_PFS ** 2 / toti_PS)) / 2 if toti_PS > 0 else -1
        if norm and tadi_norm_term:
            tadi_chi2_min = tadi_chi2_min/tadi_norm_term
        tads_chi2_min_df = pd.concat([tads_chi2_min_df,
                                      pd.DataFrame({tads_chi2_min_df.columns[0]: [comp_chr],
                                       tads_chi2_min_df.columns[1]: [tad_start],
                                       tads_chi2_min_df.columns[2]: [tad_end],
                                       tads_chi2_min_df.columns[3]: [tadi_chi2_min]})
                                     ])

        if plots_folder and (toti_PS > 0):
            logger.info(f' tadi_chi2_min={tadi_chi2_min}')
            fig.suptitle(
                f'HIC compare chromosome {comp_chr}/{comp_chr}'
                f' \n TAD: {tad_start}-{tad_end}, tot_PFS={tot_PFS}, \n tadi_chi2_min={tadi_chi2_min}')

            ax1.set_title(f'{h1_title}: \n p_i: {p_i}, \n tot_PS={tot_PS}')
            ax2.set_title(f'{h2_title}: \n p_i: {f_i}, \n tot_FS={tot_FS}')
            #plt.show()
            if not os.path.exists(plots_folder):
                try:
                    os.mkdir(plots_folder)
                except:
                    logger.warning(f'Meanwhile folder {plots_folder} was created, skip creating it!')
            fig_filename = os.path.join(plots_folder, f'{comp_filename}_{chi2_mode}_tad_{tad_start}-{tad_end}.{plot_format}')
            logger.info(f'Save figure in file: {fig_filename}')
            fig.savefig(fig_filename, format=plot_format)
            plt.close()

    tads_chi2_min_gr = pr.PyRanges(tads_chi2_min_df)
    # save tads_chi2_min_df as bigwig
    tads_chi2_min_gr.to_bed(os.path.join(plots_folder if plots_folder else '', f'{comp_filename}_tads_chi2-min_{chi2_mode}.bed'), keep=True)
    #tads_chi2_min_gr.to_bed(f'{comp_filename}_tads_chi2-min.bed', keep=True)
    #chr_sizes_gr = pr.from_dict({'Chromosome': [comp_chr], 'Start': [0], 'End': [chr_size]})
    #pr.to_bigwig(tads_chi2_min_gr, f'{comp_filename}_tads_chi2-min.bw', chr_sizes_gr)

    # final calculation
    alpha_min = (tot_PFS / tot_PS)
    chi2_min = (tot_FS - (tot_PFS ** 2 / tot_PS)) / 2

    if norm and norm_term:
        chi2_min = chi2_min / norm_term
    return chi2_min, alpha_min


def compare_hic_chromosome(hic_file, cmp_hic, hic_chrs=None, chrs=CHR_SYNONYMS, res=RESOLUTION, tads_boundary=None,
                           plots_folder=None, norm=True, chi2_mode=CHI2_MODE_LOG, hic_balance=CHI2_USE_BALANCED,
                           plot_cmap=CMAP, plot_format=PLOT_FORMAT, interaction_sites=None, lef_loading_sites=None,
                           basal_loading_factor=None,
                           lef_boundaries=None):
    """
        Compare a simulation HiC with the first HiC from a list of experimental HiCs.
        The list experimental HiCs is used to calculate the standard deviation of the average contact probability
        for a given radius.
    :param hic_chrs: list of chromosome synonyms also names of different chromosomes to be be in
    :param hic_file: hic file to be compared, if None same as chrs
    :param cmp_hic: Hi-C to compare with
    :param chrs: synonyms of chromosome to compare with
    :param res: resolution
    :param tads_boundary: TADs boundaries. Default: None, will use the whole chromosome as one TAD
    :param plots_folder: plots folder to save plots if not empty otherwise do not produce plots
    :param norm: Experimental! to normalize or not to normalize single TAD by their height (used log-lines for the statistics)
    :param chi2_mode: Chi2-min mode: linear or log
    :param hic_balance: to use ICE balanced Hi-C matrix if available, otherwise not
    """
    if hic_chrs is None:
        hic_chrs = chrs

    if os.path.exists(hic_file) and hic_file.endswith('.hdf5') and DECAY_USE_HDF5:
        hic_mat1 = get_hic(hic_file, resolution=res, balance=hic_balance)
        hic1_chr = hic_chrs[0]
        hic1_chr_size = hic_mat1.shape[0]*res
    else:
        hic1cooler, hic1_chrs = get_hic_cooler_res(hic_file, hic_chrs, res)
        # we can compare only one chromosome: we assume we have a list of chromosome synonyms
        hic1_chr = hic1_chrs[0]
        balanced = (hic1cooler.bins()['weights'] is not None) if hic_balance else hic_balance
        hic_mat1 = hic1cooler.matrix(balance=balanced).fetch(hic1_chr)
        if balanced:
            hic_mat1 = np.nan_to_num(hic_mat1)
        hic1_chr_size = hic1cooler.chromsizes[hic1_chr]

    cmp_hic_mcool = get_exp_sim_mcool(cmp_hic, chrs, res)

    hic2 = cmp_hic
    hic2cooler = cooler.Cooler(f'{cmp_hic_mcool}::/resolutions/{res}')
    exp_chr = list(set(hic2cooler.chromnames) & set(chrs))
    if exp_chr:
        exp_chr = exp_chr[0]
    else:
        logger.error(f'Experimental HiC data contains none of the chromosome names: {chrs} instead {hic2cooler.chromnames}')
        exit(1)

    # get hic file names
    h1 = hic_file[hic_file.rfind('/') + 1:]
    h2 = hic2[hic2.rfind('/') + 1:]

    # to evaluate in debug: cooler.Cooler(f'{hic_file}.{SIM_RESOLUTION}.cool').matrix(balance=False).fetch('6')
    balanced = (hic2cooler.bins()['weights'] is not None) if hic_balance else hic_balance
    hic_mat2 = hic2cooler.matrix(balance=balanced).fetch(exp_chr)  # TODO use balanced for experimental data
    if balanced:
        hic_mat2 = np.nan_to_num(hic_mat2)

    tads_pref = ('_t'+os.path.splitext(os.path.basename(tads_boundary))[0]) if tads_boundary is not None else ''
    comp_filename = f'comp_{h1}_{h2}_res{res}{tads_pref}{"_balanced" if balanced else ""}'

    chr_end = min(hic1_chr_size, hic2cooler.chromsizes[exp_chr])
    tads = read_tads_bed(tads_boundary)

    if tads is None:
        logger.info(f'Whole chromosome as a single TAD representing the whole chromosome.')
        tads = np.zeros((1, 2), dtype=int)
        tads[0, 0] = 1
        tads[0, 1] = chr_end

    (chi2_min, alpha_min) = chi2_minimization(hic_mat1, hic2cooler, chrs, res, tads, comp_filename,
                                              plots_folder=plots_folder, norm=norm, chi2_mode=chi2_mode,
                                              hic_balance=hic_balance, plot_cmap=plot_cmap, plot_format=plot_format)

    if plots_folder:
        hic_mat1_log = np.log10(hic_mat1 * alpha_min)
        hic_mat2_log = np.log10(hic_mat2)
        # split
        # fig, (ax1, ax2) = plt.subplots(1, 2)
        # fig.suptitle(f'HIC compare chromosome {hic1_chr}/{exp_chr} alpha_min:{alpha_min}')
        #
        # res_kb = res//1000
        # ax1.get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: int(x*res_kb)))
        # ax1.get_yaxis().set_major_formatter(FuncFormatter(lambda y, p: int(y*res_kb)))
        #
        # ax1.imshow(hic_mat1_log, cmap=CMAP, interpolation='nearest')  # / hic_mat1.sum()
        # ax1.set_title(h1)
        # #ax1.clim(-2.75, 0)
        #
        # ax2.imshow(hic_mat2_log, cmap=CMAP, interpolation='nearest')  # / hic_mat2.sum()
        # ax2.set_title(h2)
        #
        # ax2.autoscale(False)

        # merged
        fig = plt.figure()

        hic_merged = merge_hics(hic_mat2_log, hic_mat1_log)

        res_kb = res//1000
        plt.gca().get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: int(x*res_kb)))
        plt.gca().get_yaxis().set_major_formatter(FuncFormatter(lambda y, p: int(y*res_kb)))

        plt.imshow(hic_merged, cmap=plot_cmap, interpolation='nearest')
        plt.title(f'{h1}/{hic1_chr} vs {h2}/{exp_chr}:\n\t chi2-min={chi2_min:7.2f}, {ALPHA}={alpha_min:7.2f}')

        if DEFAULT_CLIM:
            plt.clim(DEFAULT_CLIM[0], DEFAULT_CLIM[1])  # for the demo_hic # -2.75, 0) # best found for simulations
        cbar_h = plt.colorbar()

        # hic1cooler.info['nbins'] * hic1cooler.info['bin-size']
        # chr_end = min(hic_mat1.shape[0], hic_mat2.shape[0])  # * res
        # tads = np.append(tads, chr_end)

        for tadi in tads:
            tad_start = tadi[0]/res
            tad_end = tadi[1]/res
            # ax2.plot([tad_start, tad_start+tad_end], "b-")
            logger.info(f'plot TAD {tad_start}-{tad_end}')
            # split
            # ax1.plot([tad_start, tad_end, tad_end, tad_start, tad_start],
            #          [tad_start, tad_start, tad_end, tad_end, tad_start], 'b--')
            # ax2.plot([tad_start, tad_end, tad_end, tad_start, tad_start],
            #          [tad_start, tad_start, tad_end, tad_end, tad_start], 'b--')
            # merged

            if PLOT_COMP_TADS:
                plt.plot([tad_start, tad_end, tad_end, tad_start, tad_start],
                         [tad_start, tad_start, tad_end, tad_end, tad_start], 'b--', linewidth=1, alpha=0.25)
        x_min, x_max = plt.xlim()
        # plot X-axis:boundaries; Y-axis: interaction sites, loading sites
        if interaction_sites and os.path.isfile(interaction_sites):
            int_sites_pd = pd.read_table(interaction_sites, delimiter='[,\t ;|:]', index_col=0, comment='#')
            for index, site in int_sites_pd.iterrows():
                position = site['position']//RESOLUTION
                site_len = site['length']//RESOLUTION
                plt.plot(0, position, '>', clip_on=False, color='darkgreen', markersize=1+site_len, zorder=10)

        if lef_loading_sites and os.path.isfile(lef_loading_sites):
            loading_sites_pd = pd.read_table(lef_loading_sites, delimiter='[,\t ;|:]', index_col=0, comment='#')
            for index, site in loading_sites_pd.iterrows():
                position = site['position']//RESOLUTION
                site_len = site['length']//RESOLUTION
                plt.plot(x_max, position, 'D', clip_on=False, color='cyan', markersize=1+site_len, zorder=10)
            if basal_loading_factor:
                basal_loading_factor = float(basal_loading_factor)
                if basal_loading_factor > 0.:
                    plt.vlines(x=x_max, ymin=0, ymax=x_max, color='cyan', linestyle='solid',
                               linewidth=5*basal_loading_factor)

        if lef_boundaries and os.path.isfile(lef_boundaries):
            boundaries_pd = pd.read_table(lef_boundaries, delimiter='[,\t ;|:]', index_col=0, comment='#')
            for index, site in boundaries_pd.iterrows():
                position = site['midpoint']//RESOLUTION
                plt.plot(position, 0, 'd', clip_on=False, color='slateblue', markersize=2, zorder=10)

        # plt.show()
        if not os.path.exists(plots_folder):
            try:
                os.mkdir(plots_folder)
            except:
                logger.warning(f'Meanwhile folder {plots_folder} was created, skip creating it!')
        fig_filename = os.path.join(plots_folder, f'{comp_filename}_{"balanced" if hic_balance else ""}_{chi2_mode}.{plot_format}')
        logger.info(f'Save figure in file: {fig_filename}')
        fig.savefig(fig_filename, dpi=1000, format=plot_format)
        plt.close()

    logger.info(f'chi2_minimization score for {h1} and {h2} (resolution: {res}): {chi2_min},{alpha_min} on TADs: {tads_boundary}')
    return chi2_min, alpha_min


def read_tads_bed(tads_bed_file):
    tads = None
    if tads_bed_file is not None:
        if os.path.exists(tads_bed_file):
            if pathlib.Path(tads_bed_file).suffix == '.tsv' or pathlib.Path(tads_bed_file).suffix == '.bed':
                loops_df = pr.read_bed(tads_bed_file, as_df=True)
                tads = np.zeros((loops_df.shape[0], 2), dtype=int)
                for i, loop in loops_df.iterrows():
                    tads[i, 0] = loop[1]
                    tads[i, 1] = loop[2]

        else:
            logger.info(f'Missing TADs BED format file {tads_bed_file}: single TAD mode')
    return tads


def get_exp_sim_mcool(hic, chrs, res=RESOLUTION):
    exp_sim_cool = hic if hic.endswith('.cool') else f'{hic}.cool'
    if os.path.exists(exp_sim_cool) and not os.path.exists(f'{hic}.hdf5'):  # experiment
        hic_cooler = cooler.Cooler(f'{exp_sim_cool}::/')
        root_res = hic_cooler.binsize
        factors = copy.deepcopy(RES_FACTORS)
        bisect.insort(factors, res//root_res)
        resolutions = [int(i * root_res) for i in factors]
        exp_sim_mcool = re.sub(r'.cool', f'.{resolutions[0]}.mcool', exp_sim_cool)
        try:
            res_cool = cooler.Cooler(f'{exp_sim_mcool}::/resolutions/{res}') if os.path.isfile(exp_sim_mcool) else None
        except:
            logger.info(f'Missing resolution {res} in {exp_sim_mcool} so it will be generated again')
            res_cool = None
        if not os.path.isfile(exp_sim_mcool) or not res_cool:
            exp_sim_mcool = balance_mcool(exp_sim_cool, resolutions, exp_sim_mcool)
    else:  # simulation !!!
        root_res = SIM_RESOLUTION
        factors = RES_FACTORS
        bisect.insort(factors, res//root_res)
        resolutions = [int(i * root_res) for i in factors]
        exp_sim_mcool = f'{hic}.{resolutions[0]}.mcool'
        if not os.path.exists(exp_sim_mcool):
            logger.info(f'Generating cooler files for {hic}')
            exp_sim_mcool = hic_to_mcool(hic, chrs[0], root_res, factors)
    return exp_sim_mcool


def get_hic_cooler_res(hic, chrs, res):  # , root_res=None, factors=None):
    """
    Returns a HiC cooler in a requested resolution and a list of found chromosomes out of a given list of different
    chromosomes or their synonyms.
    :param hic: input HiC file name
    :param chrs: list of different chromosomes or their synonyms
    :param res: resolution
    :return:
    """
    sim_hic_mcool = get_exp_sim_mcool(hic, chrs)  # , root_res=root_res, factors=factors)
    try:
        hic_cooler = cooler.Cooler(f'{sim_hic_mcool}::/resolutions/{res}')
    except:
        logger.warning(f'Problem to open mcool file {sim_hic_mcool} so will regenerate it!')
        os.remove(sim_hic_mcool)
        sim_hic_mcool = get_exp_sim_mcool(hic, chrs)  # , root_res=root_res, factors=factors)
        hic_cooler = cooler.Cooler(f'{sim_hic_mcool}::/resolutions/{res}')
    hic_chrs = list(set(hic_cooler.chromnames) & set(chrs))
    return hic_cooler, hic_chrs


def get_hic_cool(hic_h5, out_prefix, res=EXP_RESOLUTION):
    hic_cool = f'{out_prefix}.cool'
    org_hic_cool = f'{hic_h5}.{SIM_RESOLUTION}.cool'
    if not os.path.exists(hic_cool):
        if not os.path.exists(org_hic_cool):
            org_hic_cool = hic_to_cooler(hic_h5, chr=SIM_CHR, resolution=SIM_RESOLUTION)
        factor = res // SIM_RESOLUTION
        cmd = f'cooler coarsen -k {factor} -o {hic_cool} {org_hic_cool}'
        status, stout = subprocess.getstatusoutput(cmd)
        logger.info(f" call: {cmd}\n\t {stout}")
    return hic_cool


def plot_distance_contact_prob_decay(hic_list, hic_chrs=CHR_SYNONYMS, tads=None, exp_cool=None, output_folder=None, res=RESOLUTION,
                                     confidence=0., replace=True, chi2_mode=CHI2_MODE_LOG, format=PLOT_FORMAT):

    if hic_chrs is None:
        hic_chrs = CHR_SYNONYMS
    # plot also distance-contacts decoy plot
    if output_folder is None:  # use the folder of the first: TODO backward compatibility, could be removed later
        output_folder = os.path.dirname(hic_list[0])
    hic_names = [os.path.splitext(os.path.basename(hic))[0] for hic in hic_list]
    # we will use the simulation folder (../<radius_analyse>) of the first hic,
    # assuming they are all from the same simulation
    exp_base_name = os.path.splitext(os.path.basename(exp_cool))[0] if exp_cool else ""
    # exp_base_name = exp_base_name[0] if exp_base_name else ''
    hics = ''
    for i, hic_r in enumerate(hic_names):
        hic = hic_list[i]
        hics += '_' + hic_r + (f'.h5{"c" if DECAY_USE_HDF5_COOL else ""}'
                               if os.path.exists(hic) and hic.endswith('.hdf5') and DECAY_USE_HDF5 else '')
    tads_pref = ('_t'+os.path.splitext(os.path.basename(tads))[0]) if tads is not None else ''
    hic_decay_plot = os.path.join(output_folder,
                                  f'contact-decay_{exp_base_name}_{CHR_SYNONYMS[-1]}_vs_{hic_chrs[-1]}_{hics}{tads_pref}{f"_c{confidence:4.2f}" if confidence > 0. else ""}'
                                  f'_chi2_{chi2_mode[:3]}{"_sd" if not CHI2_USE_SEM else ""}'
                                  f'{"_zoom" if CHI2_MODE_LOG_ZOOM else ""}.{res}.{format}')
    if os.path.exists(hic_decay_plot):
        if not replace:
            logger.info(f'Distance-contact decay plot already existing: {hic_decay_plot}, so skip it')
            return hic_decay_plot
        else:
            logger.warning(f'Replacing existing Distance-contact decay plot: {hic_decay_plot}')
    logger.info(f'Distance-contact decay plot: {hic_decay_plot} ...')
    max_tad_size = chr_size // SIM_RESOLUTION
    fig = plt.figure()
    plt.axvline(CHI2_RANGE_START * SIM_RESOLUTION / res, color='gray', linewidth=1)
    plt.axvline(CHI2_RANGE_END * SIM_RESOLUTION / res, color='gray', linewidth=1)
    min_dist = CHI2_RANGE_START * SIM_RESOLUTION // res
    max_dist = CHI2_RANGE_END * SIM_RESOLUTION // res
    if chi2_mode == CHI2_MODE_LOG:
        plt.xscale('log')
        if not CHI2_MODE_LOG_ZOOM:
            min_dist = 1
            max_dist = None

    plt.yscale('log')
    lines = []
    legend = []
    cmp_hic = re.sub(r'.cool', '', exp_cool) if exp_cool else ""
    chr_i = 0
    for i, hic in enumerate(hic_list):
        if not os.path.exists(hic):
            logger.error(f'Missing file {hic}, skip it!')
            continue
        # hic_folder = os.path.dirname(hic)
        if os.path.exists(hic) and hic.endswith('.hdf5') and DECAY_USE_HDF5:
            # it's a simulation HDF5 file containing only one chromosome so let's use the first chromosome synonym
            hic_chrs_select = [hic_chrs[0]]
        else:
            hic_cooler, hic_chrs_select = get_hic_cooler_res(hic, hic_chrs, res)

        for hic_chr in hic_chrs_select:
            if os.path.exists(hic) and hic.endswith('.hdf5') and DECAY_USE_HDF5:
                (dists, probs, probs_confi_l, probs_confi_u), hic_chr_size = get_decay_distribution_hdf5(hic, resolution=res,
                                                                                                         confidence=confidence,
                                                                                                         min_dist=min_dist,
                                                                                                         max_dist=max_dist)
            else:
                dists, probs, probs_confi_l, probs_confi_u = get_decay_distribution_cool(hic_cooler, chr=hic_chr,
                                                                                         resolution=res,
                                                                                         confidence=confidence,
                                                                                         min_dist=min_dist,
                                                                                         max_dist=max_dist)
                hic_chr_size = hic_cooler.chromsizes[hic_chr]

            max_tad_size = min(max_tad_size, hic_chr_size // res)

            # calculate proper chi2_alpha
            plots_folder = os.path.join(output_folder, 'plots')
            if cmp_hic:
                chi2, alpha = compare_hic_chromosome(hic, cmp_hic, hic_chrs=[hic_chr], chrs=CHR_SYNONYMS, res=res,
                                                     tads_boundary=tads, plots_folder=plots_folder, norm=True,
                                                     chi2_mode=chi2_mode)
            else:
                chi2, alpha = 0, 1
    #        logger.info(f'COMPARE {cmp_hic} <- {hic}: {chi2}, {alpha}')
    #        chi2_rev, alpha_rev = compare_hic_chromosome(cmp_hic, hic, hic_chrs=CHR_SYNONYMS, chrs=hic_chrs, res=res, tads_csv=None,
    #                                             norm=True, chi2_mode=chi2_mode)  # , plot=True)
    #        logger.info(f'REV COMPARE {hic} <- {cmp_hic}: {chi2_rev}, {alpha_rev}')

            probs_adjust = np.array(probs) * alpha
            line, = plt.plot(dists, probs_adjust, color=DECAY_PALETTE(chr_i * DECAY_PALETTE_SHIFT), alpha=DECAY_PLOT_ALPHA)
            lines.append(line)
            if confidence != 0:  # plot confident interval otherwise not
                # see https://matplotlib.org/3.1.3/tutorials/colors/colormaps.html why I choose i*4 + 1
                plt.plot(dists, probs_confi_u, color=DECAY_PALETTE(chr_i * DECAY_PALETTE_SHIFT + 1), alpha=DECAY_PLOT_ALPHA)
                plt.plot(dists, probs_confi_l, color=DECAY_PALETTE(chr_i * DECAY_PALETTE_SHIFT + 1), alpha=DECAY_PLOT_ALPHA)
            legend.append(f'{hic_names[i]}.{hic_chr} chi2_min: {chi2:7.3f}')

            if SAVE_DECAY_PROBABILITY:  # save the exp decay probability
                f_name = f'{hic}.{res}_decay_probs_{chi2_mode[:3]}.txt'
                logger.info(f'Save simulation data from Distance-contact decay plot: {hic} to {f_name}')
                f = open(f_name, 'w')
                f.write('\n'.join([str(elem) for elem in probs]))
                f.close()

            # next chromosome
            chr_i += 1

    if exp_cool is not None:
        # plot experimental decay but correct genomic distance
        exp_cooler, exp_chr_select = get_hic_cooler_res(cmp_hic, CHR_SYNONYMS, res)
        exp_chr = exp_chr_select[0]  # supported only one chromosome
        hics_count = len(hic_list)
        # for i, exp_chr in enumerate(exp_cooler.chromnames):
        if exp_chr:
            # else:  logger.error(f'No chromosome like [{", ".join(CHR_SYNONYMS)}] was found in {exp_cooler.filename}')

            dists, probs, probs_confi_l, probs_confi_u = get_decay_distribution_cool(exp_cooler, chr=exp_chr,
                                                                                     resolution=res,
                                                                                     confidence=confidence,
                                                                                     min_dist=min_dist,
                                                                                     max_dist=max_dist)
            max_tad_size = min(max_tad_size, exp_cooler.chromsizes[exp_chr] / res)

            if SAVE_DECAY_PROBABILITY:  # save the exp decay probability
                f_name = f'{cmp_hic}.{res}_decay_probs_{chi2_mode[:3]}.txt'
                logger.info(f'Save experimental data from Distance-contact decay plot: {exp_cool} to {f_name}')
                f = open(f_name, 'w')
                f.write('\n'.join([str(elem) for elem in probs]))
                f.close()

            color = 'navy'
            # color=DECAY_PALETTE((hics_count+i) * DECAY_PALETTE_SHIFT)
            line, = plt.plot(dists, probs, color=color, alpha=DECAY_PLOT_ALPHA)
            lines.append(line)
            legend.append(f'{os.path.basename(exp_cool)}.{exp_chr}')
            if confidence != 0:  # plot confident interval otherwise not

                # see https://matplotlib.org/3.1.3/tutorials/colors/colormaps.html why I choose i*4 + 1
                plt.plot(dists, probs_confi_u, color='blue', alpha=DECAY_PLOT_ALPHA)
                plt.plot(dists, probs_confi_l, color='blue', alpha=DECAY_PLOT_ALPHA)

    # plot chi2 dist range
    logger.info(f'call get_chi2_dist_range({chi2_mode}, {res}) # calculated max_tad_size:{max_tad_size}')
    dist_range = get_chi2_dist_range(chi2_mode, res)  # TODO check if needed: , max_tad_size)
    for d in dist_range:
        plt.axvline(x=d, ymax=0.025, color='gray', linewidth=1)

    plt.legend(lines, legend)

    plt.title(f'Distance-contact decay chi2:x{chi2_mode}\n{output_folder.strip()}')
    plt.ylabel('average #contacts ~ contact probability')
    plt.xlabel(f'genomic distance in {res / 1000}kb')

    fig.savefig(hic_decay_plot, format=format)
    plt.close()
    logger.info(f' save plot: {hic_decay_plot}')
    return hic_decay_plot


def get_decay_distribution_hdf5(hic_h5, resolution, confidence=0., min_dist=1, max_dist=None):
    hic = get_hic(hic_h5, resolution)
    return get_decay_distribution(hic, confidence=confidence, min_dist=min_dist, max_dist=max_dist), hic.shape[0]


def get_hic(hic_h5, resolution=SIM_RESOLUTION, balance=CHI2_USE_BALANCED):
    if not DECAY_USE_HDF5_COOL:
        # directly from .hdf5
        with h5py.File(hic_h5, 'r') as f:
            a_group_key = list(f.keys())[0]
            logger.info(f"get_hic for {hic_h5} HDF5 Keys {f.keys()} and dimentions {a_group_key}")

            # Get the data
            data = list(f[a_group_key])
            hic = np.array(data)
            np.fill_diagonal(hic, 0)
            if resolution > SIM_RESOLUTION:  # need binning
                hic = hic_coarsen(hic, hic_h5, resolution)
    else:
        # use the root cool
        cool_file = hic_h5 + f'.{SIM_RESOLUTION}.cool'
        if not os.path.isfile(cool_file):    # and resolution == SIM_RESOLUTION:
            logger.info(f'Creating hic {cool_file}...')
            cool_file = hic_to_cooler(hic_h5, SIM_CHR, SIM_RESOLUTION)  # == cool_file
        logger.info(f'Extracting hic from {cool_file}...')
        hic_cooler = cooler.Cooler(f'{cool_file}::/')
        hic_chr = list(set(hic_cooler.chromnames) & set(SIM_CHR_SYNONYMS))[0]
        balanced = (hic_cooler.bins()['weights'] is not None) if balance else balance
        hic = hic_cooler.matrix(balance=balanced).fetch(hic_chr)
        if balanced:
            hic = np.nan_to_num(hic)  # nan -> 0
        hic = hic_coarsen(hic, hic_h5, resolution)
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


def get_decay_distribution_cool(hic_cooler, chr, resolution, confidence=0., min_dist=1, max_dist=None):
    # hic_cooler = cooler.Cooler(f'{hic_cool}::/')
    hic = hic_cooler.matrix(balance=False).fetch(chr)
    np.fill_diagonal(hic, 0)

    resolution_factor = hic_cooler.binsize / resolution
    min_dist_corrected = int(round(min_dist / resolution_factor)) if min_dist else min_dist
    max_dist_corrected = int(round(max_dist / resolution_factor)) if max_dist else max_dist
    # hic_prob = hic / np.max(hic)  # need to rethink this normalization
    dists, probs, probs_confi_l, probs_confi_u = get_decay_distribution(hic, confidence=confidence,
                                                                        min_dist=min_dist_corrected,
                                                                        max_dist=max_dist_corrected)
    dists_corrected = [int(round(d * resolution_factor)) for d in dists] if resolution_factor != 1. else dists
    # return dists_corrected, probs, probs_confi_l, probs_confi_u  # if normalized before
    probs_normed = probs  # / np.max(probs)  # wrong normalization: it has to be the same normalization factor for all
    probs_confi_l_normed = probs_confi_l  # / max(probs_confi_l) if probs_confi_l else probs_confi_l
    probs_confi_u_normed = probs_confi_u  # / max(probs_confi_u) if probs_confi_u else probs_confi_u
    return dists_corrected, probs_normed, probs_confi_l_normed, probs_confi_u_normed


def get_decay_distribution(hic, confidence=0., min_dist=1, max_dist=None):
    logger.info(f"get_decay_distribution: data.shape: ${hic.shape}, confidence:{confidence}, max_dist:{max_dist}")
    # build the cooler fields
    max_dist = min(hic.shape[0], max_dist) if max_dist else hic.shape[0]
    dists = np.arange(min_dist, max_dist)
    probs = []
    probs_confi_u = []
    probs_confi_l = []

    # hic_log = np.log10(hic)
    # plt.imshow(hic_log, interpolation='nearest', cmap='hot_r')

    for d in dists:
        avrg_prob, sd_sem_prob = average_contact_prob(hic, d)  # , plot=(d % 100 == 0), ax=plt)
        probs.append(avrg_prob)
        # calculate confidential interval
        if confidence != 0.:  # plot confident interval otherwise not
            n = max_dist - d
            test_stat = stats.t.ppf((confidence + 1) / 2, n - 1)  # t-test
            # test_stat = stats.norm.ppf((interval + 1)/2) # z-test
            if CHI2_USE_SEM:
                sem_prob = sd_sem_prob
            else:
                sem_prob = sd_sem_prob / math.sqrt(n)
            h = test_stat * sem_prob
            probs_confi_u.append(avrg_prob - h)
            probs_confi_l.append(avrg_prob + h)

    # plt.close()
    return dists, probs, probs_confi_l, probs_confi_u


def __get_num_measurements(chip_out):
    return chip_out.shape[0] if chip_out.ndim > 1 else 1


def __get_num_chain(chip_out):
    return chip_out.shape[1] if chip_out.ndim > 1 else len(chip_out)


def chip_seq(chip_out, bin_factor):
    Nmeas = __get_num_measurements(chip_out)
    Nchain = __get_num_chain(chip_out)

    bins = math.ceil(Nchain / bin_factor)
    chip_seq_binned = np.zeros(bins, dtype=float)
    for m in range(Nmeas):
        chip_seq = chip_out.iloc[m].values if isinstance(chip_out, pd.DataFrame) else chip_out[
            m] if chip_out.ndim > 1 else chip_out
        for i in range(bins):
            chip_seq_binned[i] += np.sum(chip_seq[(i * bin_factor):min((i + 1) * bin_factor, Nchain)])
    chip_seq_binned = chip_seq_binned / Nmeas

    return chip_seq_binned


def chip_out_to_bedgraph(chip_out_file: str, bed_graph_file: str = None, chrom=SIM_CHR, resolution: int = SIM_RESOLUTION):
    if not bed_graph_file:  # default
        bed_graph_file = f'{os.path.splitext(chip_out_file)[0]}.bedGraph'

    bin_factor = resolution // SIM_RESOLUTION

    chip_out = pd.read_csv(chip_out_file, delim_whitespace=True, encoding='utf-8')
    chip_seq_binned = chip_seq(chip_out, bin_factor=bin_factor)  # bin_factor = 5 : 10kb = 5*2kb
    nbins = __get_num_chain(chip_seq_binned)

    chip_seq_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'value'], index=None)
    chip_seq_df['start'] = np.arange(0, nbins*resolution, resolution)
    chip_seq_df['end'] = np.arange(resolution, (nbins+1)*resolution, resolution)
    chip_seq_df['value'] = chip_seq_binned
    chip_seq_df['chrom'] = chrom

    exp_chip_pd = chip_seq_df.to_csv(bed_graph_file, sep='\t', header=False, index=False)
    return exp_chip_pd


def read_boundary_pos(boundary_file, Nchain):
    boundary_pos = np.zeros(Nchain, dtype=float)
    if boundary_file is not None:
        if os.path.exists(boundary_file):
            boundary_b_positions = pd.read_csv(boundary_file, header=0)[['b-position']]['b-position']
            for i in boundary_b_positions:
                boundary_pos[i] = 1.
        else:
            logger.error(f'Boundary file not found: {boundary_file}')
    return boundary_pos


def plot_chip_seq(chip_out_file, exp_chip, boundary, resolution=SIM_RESOLUTION, correlation=DEFAULT_CHIP_CORRELATION, plot=True,
                  replace=True):
    """
    Returns correlation coefficient and optionally produce a comparison plot
    :param chip_out_file: the original simulation chip-seq.out file
    :param exp_chip: the experimental Chip-seq file in the same resolution as the simulation (2kb)
    :param boundary: the boundary file used during the simulation
    :param resolution: resolution to downscale the Chip-seq
    :param correlation: the comparative correlation to be used
    :param plot: to plot or not in a file
    :param replace: to replace existing plot file
    :return: the correlation as: (correlation coefficient, p-value, L2 distance (simulation - experiment)).
    """
    bin_factor = SIM_RESOLUTION // resolution
    chip_out = pd.read_csv(chip_out_file, delim_whitespace=True, encoding='utf-8')
    output_folder = os.path.dirname(chip_out_file)
    measurements = chip_out.shape[0]  # + 1  # +1: for the first that has no chip
    # chip_seq_m = chip_out.iloc[-1]  # take the last measurement
    sim_chip_seq = chip_seq(chip_out, bin_factor=bin_factor)  # bin_size = 5 : 10kb = 5*2kb
    sim_chip_seq_normed = sim_chip_seq / sum(sim_chip_seq)
    # read BedGraph file format
    if exp_chip:
        exp_chip_pd = pd.read_csv(exp_chip, delim_whitespace=True, names=['chrom', 'start', 'end', 'value'],
                                  encoding='utf-8')
        exp_chip_x_pd = exp_chip_pd[exp_chip_pd.chrom.isin(CHR_X_SYNONYMS)]
        exp_chip_out = exp_chip_x_pd[['value']].values.transpose()
        exp_chip_x = chip_seq(exp_chip_out, bin_factor=bin_factor)
        # sum_x = sum(exp_chip_x['value'])
        exp_chip_x_normed = exp_chip_x / sum(exp_chip_x)
    boundary_pos = read_boundary_pos(boundary, Nchain=chip_out.shape[1])
    boundary_chip = chip_seq(boundary_pos, bin_factor=bin_factor)
    boundary_chip_normed = boundary_chip * len(np.where(boundary_chip > 0.)[0]) \
        / (max(exp_chip_x) * sum(boundary_chip) * 100) if exp_chip else boundary_chip
    if exp_chip:
        chip_l2 = np.linalg.norm(exp_chip_x_normed - sim_chip_seq_normed)
        if correlation == 'spearmanr':
            corr = stats.spearmanr(exp_chip_x_normed, sim_chip_seq_normed)
            chip_correlation = (corr.correlation, corr.pvalue, chip_l2)
        elif correlation == 'pearsonr':
            corr = stats.pearsonr(exp_chip_x_normed, sim_chip_seq_normed)
            chip_correlation = (corr[0], corr[1], chip_l2)
    else:
        chip_correlation = 0

    if plot:
        size_kb = resolution // 1000
        plot_file_name = os.path.join(output_folder, f'chip-seq_{measurements}_{size_kb}kb.png')
        if os.path.exists(plot_file_name):
            if not replace:
                logger.info(f'Chip-seq plot already existing: {plot_file_name}, so skip it')
                return chip_correlation
            else:
                logger.warning(f'Replacing existing Chip-seq plot: {plot_file_name}')
        fig = plt.figure()
        if exp_chip:
            exp_chip_line = plt.plot(exp_chip_x_normed, color='b', linewidth=1, alpha=DECAY_PLOT_ALPHA)
            chip_correlation_title = f'{correlation}: {chip_correlation[0]:8.5f}, p-val:{chip_correlation[1]:8.5f}, L2:{chip_correlation[2]:8.5f}'
        else:
            chip_correlation_title = ""
        sim_chip_line = plt.plot(sim_chip_seq_normed, color='r', linewidth=1, alpha=DECAY_PLOT_ALPHA)
        boundary_chip_line = plt.plot(boundary_chip_normed, color='black', linewidth=1, alpha=DECAY_PLOT_ALPHA)
        sim_folder = output_folder   # .split('/')[-2]
        plt.title(f'ChIP-seq for simulation \n {sim_folder} snapshots: {measurements} \n{chip_correlation_title}')
        plt.ylabel('average #bound LEFs')
        plt.xlabel(f'genomic position in {size_kb}kb')
        plt.legend(('experimental', 'simulation', 'boundary'))
        fig.savefig(plot_file_name)
        plt.close()

    return chip_correlation


def merge_hics(hic_mat_upper, hic_mat_lower):
    if len(hic_mat_upper) != len(hic_mat_lower):
        logger.warning(f'Comparing different size of Hi-C matrices: '
                       f'size(hic_mat_upper)={len(hic_mat_upper)} != size(hic_mat_upper)={len(hic_mat_lower)}.'
                       f'It will be used the smaller size for the merged plot!')
    for index_row in range(min(len(hic_mat_upper), len(hic_mat_lower))):
        hic_mat_upper[index_row][:index_row] = hic_mat_lower[index_row][:index_row]

    np.fill_diagonal(hic_mat_upper, 0)

    return hic_mat_upper
