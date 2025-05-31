#! /usr/bin/env python

import argparse
import configparser
import copy
import fnmatch
import logging
import os
import threading
import re

import numpy as np
import pandas as pd
import sys
from importlib.resources import files, as_file

from py3dpolys_le import hic_analysis as ha
from py3dpolys_le.job_runner import CfgJobRunner, CFG_SECTION_3DPOLYS_LE

EXP_COOL_AS_STATS = '.'

# Initialization
logger = logging.getLogger(__name__)


def cli_parser():
    p = argparse.ArgumentParser(description=f'''Running 3dpolys_le_runner.
    Run one of the following batch commands:
    
    •	grid_nlef_km: Perform a series of simulations rotating over a range of Nlef (--nlef_list) and 
    km (--km_list) parameter values. In combination with a --replace parameter.
    
    •	new_stats: Perform comparative statistical analysis on all entries, representing simulations,  in a given 
    sim_stats.tsv file (--stats_file) comparing simulations with the already used or a new given (--exp_cool) experimental 
    data. In combination with a –-replace parameter.
    
    •	decay_plots: Produce distance-contact-decay plots on all entries in a given comparative statistics file (--stats_file) 
    and compare simulations with the already used or a new given (--exp_cool) experimental data. In combination with  
    --resolution, –-replace and --threading parameters.
    
    •	chip_seq_plots: Produce ChIP-seq plots on all entries in a given sim_stats.tsv file (--stats_file) and compare 
    simulations' ChIP-seq with the already used or a new given (--exp_chip) experimental data. In combination with  
    –-replace and --threading parameters.
    
    •	contact_radius_analysis: Perform contact analysis step (Hi-C data extraction) for new given parameters on all 
    entries in a given sim_stats.tsv file (--stats_file). In combination with --replace and --threading parameters.
    
    •	multi_decay_plot: Produce a multi distance-contact-decay plot for a single simulation’s output (--output_folder) 
    including all its subfolders (for different contact-radius analysis) and compare them with a given (--exp_cool) 
    experimental data. In combination with the –resolution and –-replace parameters.
    
    •	multi_decay_exps_plot: Produce multi distance-contact-decay plot for a given list of datasets (----exp_cools, 
    experimental, and simulations) and compare one of their chromosomes (--cmp_chrs) with a chromosome (--hic_chrs) from the 
    first dataset in the list. In combination with the –resolution and –-replace parameters.
    
    •	run: Run a single simulation including extracting Hi-C matrices (‘analyse’ step), collection all comparative 
    statistics (sim_stats.tsv file), and distance-contact-decay plots for all measurements.
    
    ''', epilog='''DISCLAIMER: 
    As almost not sufficient tests prove the correctness of all possible parameter combinations (no comprehensive test coverage), 
    please check your output data and log files, and make sure all went as you have expected.
    ''', formatter_class=argparse.RawDescriptionHelpFormatter)  # RawTextHelpFormatter
    p.add_argument("run_command", help="Run batch command.",
                   choices=['grid_nlef_km', 'new_stats', 'decay_plots', 'chip_seq_plots', 'contact_radius_analysis',
                            'multi_decay_plot', 'multi_decay_exps_plot', 'run'])
    p.add_argument("-bd", "--boundary_direction", default=0, type=int,
                   help="Impermeability direction applied to all boundaries: -1:opposite direction, 0:both, 1:same direction."
                        " Default: 0")
    # TODO add -bf and look for others missing; remove unneeded and confusing
    p.add_argument("-b", "--boundary", default="",
                   help="Boundary sites file in CSV format used in a simulation.")
    p.add_argument("-lls", "---lef_loading_sites", default="",
                   help="<loop extrusion loading sites file> LEFs loading sites file in a csv format with the following "
                        "columns: name,position,length,factor. Default: if not given, the whole polymer.")
    p.add_argument("--hic_chrs", nargs='*', default=None,
                   help="Synonyms of the chromosome from Hi-C matrixes to be compared. "
                        "Default X, set by constant hic_analysis.py::CHR_SYNONYMS = CHR_X_SYNONYMS.")
    p.add_argument("--cmp_chrs", nargs='*', default=None,
                   help="Synonyms of the Hi-C chromosome with which simulation (or other data) to be compared. "
                        "Overwrites the default constant hic_analysis.py::CHR_SYNONYMS!")
    p.add_argument("-z", "--z_loop", action='store_true', help="Allow z_loop for LEFs move in a simulation.")
    p.add_argument("-u", "--unidirectional", action='store_true',
                   help="Unidirectional mode for LEFs move otherwise bidirectional.")
    p.add_argument("-im", "--init_mode", default='',
                   help="Initial folding mode: h for helices like, z for zigzag like polymer state"
                        ", s=<sim_out_folder> to continue from a finished simulation output folder. "
                        "Default: z.")
    p.add_argument("-t", "--tads_boundary", default="",
                   help="TADs boundary file used to calculate the chi-2-min score in the same format as the "
                        "boundary sites file used for simulation. Default: no boundaries equivalent to"
                        "the whole chromosome seen as a single TAD.")
    p.add_argument("-e", "--exp_cool", default="",
                   help=f"Experimental cooler (.cool) file with which a simulation data to be compared. "
                        f"When used in combination with a new_stats command, if the value is '{EXP_COOL_AS_STATS}', "
                        f"it will use values stored in a given sim_stats.tsv file (--stats_file).")
    p.add_argument("-es", "--exp_cools", nargs='+',   # TODO maybe rename it as it is used for mixed list of datasets
                   help=f"In combination with a 'multi_decay_exps_plot' command sets the list of experimental cooler files "
                        f"(.cool) or simulation Hi-C HDF5 files, and calculates chi2-min scores to compare them all with "
                        f"the first one in the list in a multi distance-contact-decay plot.")
    p.add_argument("-ec", "--exp_chip",
                   help=f"Experimental ChIP-seq file (in .bedGraph format) with which simulation’s ChIP-seq file "
                        f"({ha.CHIP_OUT}) to be compared.")
    p.add_argument("-res", "--resolution", default=ha.RESOLUTION, type=int,
                   help="Resolution for distance-contact-decay plots, and to downscale chip-seq plots.")
    p.add_argument("-i", "--input_cfg", default="./input.cfg", help="Input.cfg file used from a simulation.")
    p.add_argument("-l", "--nlef", help="Nlef value to use in a simulation.", type=int)
    p.add_argument("-m", "--km", help="km value to use in a simulation.", type=float)
    # 0.00054, 0.00162, 0.0027, 0.00378, 0.00486, 0.00648
    # 10,      30,      50,     70,      90,      120 kb/min
    # 5.4e-4 = 10kb/min = 167bp/s

    p.add_argument("-f", "--stats_file", default="./sim_stats.csv", help="Simulation statistics' repository file.")
    p.add_argument("--stats", action='store_true',
                   help="In combination with a new_stats command to run statistical analysis (3dpolys_le_stats.py) only for "
                        "entries in a given sim_stats.csv file (--stats_file) missing statistics plots.")
    p.add_argument("--all_stats", action='store_true',
                   help="In combination with a new_stats command to run statistical analysis (3dpolys_le_stats.py) for all "
                        "entries in a given sim_stats.csv file (--stats_file).")
    # TODO unify/clarify what to use better: only -r, or -lr
    p.add_argument("-r", "--radius_contact", default=0., type=float,
                   help="Contact radius in lattice units (1=70nm) to run a single 'analyse' step "
                        "(py3dpolys_le program module) for extracting Hi-C matrixes.")
    p.add_argument("-lr", "--list_contact_radii", nargs='+', default=['2.84'],  # ['1.42', '2.13', '2.84', '3.55', '4.26'],  #, '5.0p'
                   help="list of contact radii to be used by a 'contact_radius_analysis' batch command, "
                        "a <p> at the end denote use of contact radius probability mode.")
    p.add_argument("--replace", help="Replace of any existing output data files (from the previous run) otherwise skip the "
                                     "simulation and continue with the next.", action='store_true')
    p.add_argument("--threading", help="Use Python multi-threading. In combination with a decay_plots batch command.",
                   action='store_true')
    # p.add_argument("--simultaneously", default=1, help="Run simultaneously <n> simulation jobs.", type=int)
    p.add_argument("-cp", "--contact_probability", action='store_true',
                   help="Together with the monomer contact radius use contact radius probability: (1 - r^2 / max_r^2), "
                        "where <max_r> is the value of the radius_contact parameter.")
    p.add_argument("--correlation", default='spearman', choices=['spearmanr', 'pearsonr'],
                   help="Correlation method to use to compare Chip-seq profiles.")
    p.add_argument("-o", "--output_folder", default="", help="Simulation output folder. If empty, it will be autogenerated.")
    p.add_argument("--cmd_run_file", default="", help="File where to save all commands.")
    p.add_argument("-a", "--analysis_folder", default="", help="Analysis output folder: Hi-C, Chip-Seq in silico.")
    p.add_argument("--nlef_list", nargs='+', default=[50, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200],
                   help="List of Nlef values. In combination with grid_nlef_km and grid_nlef_km_dirlef commands.", type=int)
    p.add_argument("--km_list", nargs='+', default=[5.4e-4, 3*5.4e-4, 5*5.4e-4, 7*5.4e-4, 9*5.4e-4, 12*5.4e-4],
                   # 0.00054, 0.00162, 0.0027, 0.00378, 0.00486, 0.00648
                   # 10,      30,      50,     70,      90,      120 kb/min
                   # 5.4e-4 = 10kb/min = 167bp/s
                   help="List of km values. In combination with grid_nlef_km commands.", type=float)
    # TODO refactor help to group params with commands
    # TODO load params properly for new_stats: CLI parameters over in file saved values
    # p.add_argument("--no_overwrite", help="Overwrite old files", action='store_false')
    return p


class DccExtrusionArgs:
    input_cfg: str
    output_folder: str
    analyse: str
    boundary: str
    lef_loading_sites: str
    tads_boundary: str
    stats_file: str
    exp_cool: str
    exp_chip: str
    nlef: int
    km: float
    radius_contact: float
    contact_probability: bool
    boundary_direction: int
    z_loop: bool
    unidirectional: bool
    init_mode: str
    stats: bool
    all_stats: bool
    # simultaneously: int
    cmp_chrs: list
    resolution: int
    _config = configparser.ConfigParser()

    def __init__(self, stats=False, all_stats=False, boundary="", lef_loading_sites="", input_cfg="./input.cfg",
                 tads_boundary="",
                 stats_file="./py3dpolys_le_stats.csv",
                 exp_cool="",
                 exp_chip="",
                 nlef=0, km=0, radius_contact=0, contact_probability=None,
                 boundary_direction=None, z_loop=None, unidirectional=None,
                 init_mode='',
                 output_folder='', analyse='', cmp_chrs=None, resolution=ha.RESOLUTION):  # , simultaneously=1
        # input arguments values overwriting configuration values (input.cfg)
        self.input_cfg = input_cfg
        if os.path.exists(self.input_cfg):
            self._config.read(self.input_cfg)
        else:
            logger.warning(f'Input configuration file {self.input_cfg} does not exists or not accessible, '
                           f'all parameters will be read from the CLI!')

        self.output_folder = output_folder
        self.analyse = analyse
        self.boundary = boundary if boundary or not os.path.exists(self.input_cfg) else self.get_property('boundary')
        if self.boundary and not os.path.exists(self.boundary):
            logger.error(f'Boundary file {self.boundary} not found!')
            exit(1)
        self.lef_loading_sites = lef_loading_sites if lef_loading_sites or not os.path.exists(self.input_cfg) else self.get_property('lef_loading_sites')
        if self.lef_loading_sites and not os.path.exists(self.lef_loading_sites):
            logger.error(f'LEF loading-sites file {self.lef_loading_sites} not found!')
            exit(1)
        self.tads_boundary = tads_boundary if tads_boundary or not os.path.exists(self.input_cfg) else self.get_property('tads_boundary')
        if self.tads_boundary and not os.path.exists(self.tads_boundary):
            logger.error(f'TADs-boundary file {self.tads_boundary} not found!')
            exit(1)
        self.stats_file = stats_file
        self.exp_cool = exp_cool if exp_cool or not os.path.exists(self.input_cfg) else self.get_property('exp_cool')
        self.exp_chip = exp_chip  # TODO probably remove
        # self.exp_ins_score = exp_ins_score
        self.nlef = nlef if nlef or not os.path.exists(self.input_cfg) else int(self.get_property('Nlef'))
        self.km = km if km or not os.path.exists(self.input_cfg) else float(self.get_property('km'))
        self.radius_contact = radius_contact if radius_contact or not os.path.exists(self.input_cfg) else float(self.get_property('radius_contact'))
        self.contact_probability = contact_probability if contact_probability is not None or not os.path.exists(self.input_cfg) else self.get_property('contact_probability', False)  # TODO probably remove
        self.boundary_direction = boundary_direction if boundary_direction is not None or not os.path.exists(self.input_cfg) else int(self.get_property('boundary_direction'))
        self.z_loop = z_loop if z_loop or not os.path.exists(self.input_cfg) or not os.path.exists(self.input_cfg) else self.get_property('z_loop').lower() in ['true', '1', 't', 'y', 'yes']
        self.unidirectional = unidirectional if unidirectional is not None or not os.path.exists(self.input_cfg) else self.get_property('unidirectional').lower() in ['true', '1', 't', 'y', 'yes']
        self.init_mode = init_mode if init_mode or not os.path.exists(self.input_cfg) else self.get_property('init_mode')
        self.stats = stats
        self.all_stats = all_stats
        self.cmp_chrs = cmp_chrs if cmp_chrs or not os.path.exists(self.input_cfg) else re.split('\\s*;\\s*|\\s*,\\s*|\\s+', self.get_property('cmp_chrs'))
        self.resolution = resolution
        if not self.cmp_chrs:
            # if nothing else final default: c.elegans chrX
            self.cmp_chrs = ha.CHR_SYNONYMS
        else:
            ha.CHR_SYNONYMS = self.cmp_chrs
        # self.simultaneously = simultaneously
        if self.exp_cool:
            if not os.path.exists(self.exp_cool):
                logger.error(f'Experimental cool file {self.exp_cool} not found!')
                exit(1)
            else:
                # make sure needed .mcool are created to avoid concurrency problem
                # i.e. while new_stats with new exp_cool
                logger.info(f'Ensure .mcool file of {self.exp_cool} is created...')
                cmp_hic = re.sub(r'.cool', '', self.exp_cool)
                ha.get_exp_sim_mcool(cmp_hic, self.cmp_chrs, self.resolution)

    def get_property(self, name, default=''):
        try:
            value = self._config.get(CFG_SECTION_3DPOLYS_LE, name)  # .strip()
        except (configparser.NoOptionError, configparser.NoSectionError) as e:
            value = default
        return value

    def default_analysis_subfolder(self) -> str:
        cp = 'p' if self.contact_probability else ''
        return f'r{self.radius_contact:4.2f}{cp}' if self.radius_contact > 0 else ''

    def default_analysis_folder(self) -> str:
        return os.path.join(self.output_folder, self.default_analysis_subfolder())

    def default_output_folder(self) -> str:
        """
        Constructs a default output folder based on the running parameters
        also add a sub-folder for the specified contact radius.

        :return: default output folder
        :rtype: str
        """
        z_opt = '_z-loop' if self.z_loop else ''
        u_opt = '_unidir' if self.unidirectional else ''
        im_opt = f'_im-{self.init_mode[0:1]}' if self.init_mode else ''
        return os.path.join('', f'out-Nlef{self.nlef}-km{self.km:g}-bd{self.boundary_direction}{im_opt}{z_opt}{u_opt}')


def is_analysis_output_folder(folder: str):
    chip_out_file = os.path.join(folder, ha.CHIP_OUT)
    return os.path.exists(folder) and os.path.isdir(folder) \
           and os.path.join(folder, ha.CHIP_OUT) \
           and os.path.join(folder, ha.XYZCONFIG_OUT)


def is_simulation_output_folder(folder: str) -> bool:
    return os.path.exists(folder) and os.path.isdir(folder) and os.listdir(folder) \
           and os.path.exists(os.path.join(folder, ha.CONFIG_OUT)) \
           and os.path.exists(os.path.join(folder, ha.CONTACT_OUT)) \
           and os.path.exists(os.path.join(folder, ha.NLEF_OUT)) \
           and os.path.exists(os.path.join(folder, ha.PROCESS_OUT)) \
           and os.path.exists(os.path.join(folder, ha.DR_OUT))


def str2bool(b):
    # try:
    return b.lower() in ('true', 't', 'yes', 'y', "1")
    # except:  # should not be needed but in some cases to recover from corrupted file
    #    logger.error(f'Failed to convert {b} to bool!')
    #    return False


class DccExtrusionRunner:
    _running_jobids = []
    _job_runner: CfgJobRunner

    def __init__(self, job_runner: CfgJobRunner):
        self._job_runner = job_runner

    @staticmethod
    def read_stats_file(stats_file):
        # read float as string to avoid rounding errors if decide to save it back
        return pd.read_csv(stats_file, delimiter=',', encoding='utf-8', header=0,
                           dtype={'boundary': str, 'boundary_direction': str, 'km': str, 'radius_contact': str,
                                  'chi2_log': str, 'alpha_log': str, 'chi2_lin': str, 'alpha_lin': str}).replace(np.nan, '', regex=True)
        # , engine='python')

    def run(self, dcc_args: DccExtrusionArgs, stats_only=False, radii=[], dep_jobid: str = None,
            replace: bool = False) -> str:
        '''
        Runs the main py3dpolys_le with given parameters either to run simulation or analysis afterwards.
        :param dcc_args: simulation parameters
        :param stats_only: do stats on analysis # TODO revise where to keep stats_only: use dcc_args.stats or not
        :param radii: list of contact radii to perform contact radius analysis after the simulation is done.
                    If the last character is not a digit it will be used the contact radius probability
                    formula: (1 - r^2 / max_r^2)
        :param dep_jobid: if present wait for this job to finish, probably the simulation for the analysis
                            otherwise wait for the previous running worker job probably another simulation.
                            If present but whitespace run simultaneously.
        :param replace: to replace or not an existing simulation or analysis
        :return: the job id if worker job was run
        '''
        prev_jobid = dep_jobid if dep_jobid else self._running_jobids[-1] if len(self._running_jobids) > 0 else ""
        cp = '-cp' if dcc_args.contact_probability else ''
        z_loop = '-z' if dcc_args.z_loop else ''
        u_opt = '-u' if dcc_args.unidirectional else ''
        init_mode = f'-im:{dcc_args.init_mode}' if dcc_args.init_mode else ''
        r = dcc_args.radius_contact
        # logger.info(f'r:{r}, cp:{cp}')
        r_opt = f'-r:{r} {cp}' if r > 0 else ''
        r_opt_py = f'-r {r} {cp}' if r > 0 else ''

        boundary = dcc_args.boundary
        input_cfg = dcc_args.input_cfg
        stats_file = dcc_args.stats_file
        exp_cool = dcc_args.exp_cool
        nlef = dcc_args.nlef
        km = dcc_args.km
        a_opt = ''
        if dcc_args.analyse:
            a_opt = f'-a:{dcc_args.analyse}'

        if not dcc_args.output_folder:  # make sure you have one
            dcc_args.output_folder = dcc_args.default_output_folder()

        jobid = ''
        if (not dcc_args.stats) and (not dcc_args.all_stats) and (not stats_only):
            # TODO check the usage of dcc_args.stats and  dcc_args.all_stats and clean up
            # check if already exist and proceed if allowed to replace
            run_sim_or_analysis = True
            if not dcc_args.analyse:
                if is_simulation_output_folder(dcc_args.output_folder):
                    if replace:
                        logger.warning(f'Simulation output folder {dcc_args.output_folder} '
                                       f'already exists and data files '
                                       f'({ha.CONFIG_OUT},{ha.CONTACT_OUT},{ha.DR_OUT},{ha.NLEF_OUT},{ha.PROCESS_OUT})'
                                       f' are there and will be replaced!')
                    else:
                        logger.warning(
                            f'Simulation output folder {dcc_args.output_folder} already exists and data files'
                            f'({ha.CONFIG_OUT},{ha.CONTACT_OUT},{ha.DR_OUT},{ha.NLEF_OUT},{ha.PROCESS_OUT})'
                            f' are there so simulation will be skipped!')
                        run_sim_or_analysis = False

            if dcc_args.analyse:
                if is_analysis_output_folder(dcc_args.analyse):
                    if replace:
                        logger.warning(f'Analysis output folder {dcc_args.analyse} exists and data files '
                                       f'({ha.CHIP_OUT}, {ha.XYZCONFIG_OUT}) are there and will be replaced!')
                    else:
                        logger.warning(f'Analysis output folder {dcc_args.analyse} exists and data files  '
                                       f'({ha.CHIP_OUT}, {ha.XYZCONFIG_OUT}) are there so analysis will be skipped!')
                        run_sim_or_analysis = False

            if run_sim_or_analysis:
                # workaround for Slurm: sbatch expect shell script. For other could be unneeded
                # not nice but maybe can be improved later
                cmd_prefix = self.get_cmd_prefix(dcc_args, mpirun=True)
                boundary_opt_f = f'-b:{boundary}' if boundary else ''   # fortran
                lef_loading_sites_opt_f = f'-lls:{dcc_args.lef_loading_sites}' if dcc_args.lef_loading_sites else ''
                boundary_direction_opt_f = f'-bd:{dcc_args.boundary_direction}' if dcc_args.boundary_direction is not None else ''
                cmd = f"{cmd_prefix} 3dpolys_le " \
                      f"-o:{dcc_args.output_folder} --km:{km} --nlef:{nlef} " \
                      f"{boundary_opt_f} {lef_loading_sites_opt_f} " \
                      f"{boundary_direction_opt_f} " \
                      f"{z_loop} {u_opt} {init_mode} {a_opt} {r_opt} {input_cfg}"
                if dcc_args.analyse:
                    jobid = self._job_runner.run_cmd(cmd, prev_jobid, profile='analysis')
                else:
                    jobid = self._job_runner.run_cmd(cmd, prev_jobid, profile='sim')

        if not dep_jobid and jobid:
            self._running_jobids.append(jobid)

        # ######## contact_radius_analysis ########
        stats_dep_jobid = jobid
        if not dcc_args.analyse and not stats_only:
            dcc_args_analysis = copy.deepcopy(dcc_args)
            stats_dep_jobid = self.sim_contact_radius_analysis(dcc_args_analysis, radii=radii, dep_jobid=jobid)

        # ########  Statistics ########
        #  calc stats and store them if not done already: Chip.out exists but .mcool doesn't exists
        # if dcc_args.analyse and not stats_only:
        analyse_folder = dcc_args.analyse

        hic_mcool = fnmatch.filter(os.listdir(analyse_folder), 'hic_*.mcool') if os.path.exists(analyse_folder) else []
        if analyse_folder or (
                os.path.exists(os.path.join(analyse_folder, ha.CHIP_OUT)) &
                ((len(hic_mcool) == 0) or dcc_args.all_stats or stats_only)):
            s_cmp_chrs = f'--cmp_chrs {" ".join(dcc_args.cmp_chrs)}' if dcc_args.cmp_chrs is not None else ''

            t_opt = f"-t {dcc_args.tads_boundary}" if dcc_args.tads_boundary else ""
            boundary_opt_py = f'-b {boundary}' if boundary else ''  # python
            boundary_direction_opt_py = f'-bd {dcc_args.boundary_direction}' if dcc_args.boundary_direction is not None else ''
        # for LOCAL use something like : #
            cmd_prefix = self.get_cmd_prefix(dcc_args, mpirun=False)
            cmd = f"{cmd_prefix} 3dpolys_le_stats " \
                  f"-o {dcc_args.output_folder} -a {analyse_folder} --km {km} --nlef {nlef} -e {exp_cool} " \
                  f"{boundary_opt_py} {boundary_direction_opt_py} " \
                  f"{t_opt} {r_opt_py} " \
                  f"-i {input_cfg} -f {stats_file} {s_cmp_chrs} "
                #  f"-eis {dcc_args.exp_ins_score} " \
            jobid = self._job_runner.run_cmd(cmd, stats_dep_jobid, profile='stats')  # need to wait for before continue with multi-decay

        return jobid

    def get_cmd_prefix(self, dcc_args, mpirun=False):
        # TODO find a better way to separate the Slurm problem
        # cmd_sh = resources.path('py3dpolys_le.bin', 'cmd.sh')  # TODO try 'cmd.sh'(setup.py), maybe not working
        with as_file(files('py3dpolys_le.bin') / 'cmd.sh') as cmd_sh_path:
            cmd_sh = str(cmd_sh_path)
            container_prefix = self._job_runner.get_property(profile='', name='container_prefix')
        if container_prefix:
            cmd_sh = os.path.join(dcc_args.output_folder, 'cmd.sh')
            if not os.path.exists(cmd_sh):
                if not os.path.exists(dcc_args.output_folder):
                    os.mkdir(dcc_args.output_folder)
                with open(cmd_sh, 'w') as f:
                    f.write('#! /bin/bash\n')
                    f.write('"$@"\n')
        cmd_prefix = f"{cmd_sh}{' mpirun' if mpirun else ''} {container_prefix}"  # TODO make 'mpirun ' parameterized
        return cmd_prefix

    def analysis_stats(self, dcc_args: DccExtrusionArgs, new_stats=False, exp_cool=None):
        """ deprecated for new_stats=False, replaced by contact_radius_analysis()
        Redo the analysis step including the stats (sim_stat.csv) or only regenerating the stats
        :param exp_cool: the experimental HiC cool file to compare with.
                            If none it will be take as it is in the sim_stats.csv file
        :param dcc_args: the running py3dpolys_le and analysis  arguments
        :param new_stats: whether to regenerate all statistics again (sim_stat.csv) and keep a backup copy of the old one
        :return: updated or new sim_stats.csv file and plots if requested
        """
        if new_stats:
            stats_file_bak = dcc_args.stats_file + '.bak'
            bak_i = 0
            while os.path.exists(stats_file_bak):
                stats_file_bak = dcc_args.stats_file + f'.bak{bak_i}'
                bak_i += 1
            os.rename(dcc_args.stats_file, stats_file_bak)
            sim_stats_pd = DccExtrusionRunner.read_stats_file(stats_file_bak)
            # a new args.stats_file will be generated
        else:
            sim_stats_pd = DccExtrusionRunner.read_stats_file(dcc_args.stats_file)

        columns = list(sim_stats_pd)
        done_sim_out_folders = []
        # for index, sim in enumerate(sim_stats_pd.itertuples(index=False)):   # sim_stats_pd.iterrows():
        try:
            for index, sim in sim_stats_pd.iterrows():
                logger.info(
                    f"index:{index}, sim: {sim['sim_out_folder']}, {sim['sim_hic_file']}, {sim['nlef']}, {sim['km']}")
                dcc_args_analysis = copy.deepcopy(dcc_args)
                if exp_cool and exp_cool != EXP_COOL_AS_STATS:
                    dcc_args_analysis.exp_cool = exp_cool
                else:
                    dcc_args_analysis.exp_cool = sim['exp_cool']
                # dcc_args_analysis.exp_chip = sim['exp_chip']  # TODO handle exp_chip properly
                dcc_args_analysis.boundary = sim['boundary']
                dcc_args_analysis.boundary_direction = sim['boundary_direction'] if columns.__contains__(
                    'boundary_direction') else 0
                dcc_args_analysis.nlef = sim['nlef']
                dcc_args_analysis.km = sim['km']
                rcp = sim['radius_contact']   # in the sim_stat.csv is saved in this format TODO use proper column name
                cp, r = self.split_radius_contact_probability(rcp)
                dcc_args_analysis.radius_contact = r
                dcc_args_analysis.contact_probability = cp
                dcc_args_analysis.output_folder = sim['sim_out_folder']
                hic_h5 = sim['sim_hic_file']
                dcc_args_analysis.analyse = os.path.dirname(hic_h5)

                if not is_simulation_output_folder(dcc_args_analysis.output_folder):
                    logger.error(
                        f'Output folder {dcc_args_analysis.output_folder} not existing or simulation data files '
                        f'({ha.CONFIG_OUT}, {ha.CONTACT_OUT}, {ha.DR_OUT}, {ha.NLEF_OUT}, {ha.PROCESS_OUT}) are missing so skip it!')
                    continue

                # stats_only has higher priority and will be checked first
                jobid = self.run(dcc_args_analysis, stats_only=new_stats)

                # run multi-decay-plot once for the sim_output_folder
                if new_stats and (dcc_args_analysis.output_folder not in done_sim_out_folders):
                    self.run_multi_decay_plot(dcc_args_analysis, dep_jobid=jobid)
                    done_sim_out_folders.append(dcc_args_analysis.output_folder)
        except:
            if new_stats and not os.path.exists(dcc_args.stats_file):
                # roll back from backup
                logger.error(f' Failed to process file {dcc_args.stats_file}: {sys.exc_info()[0]}')
                os.rename(stats_file_bak, dcc_args.stats_file)
                raise

    def contact_radius_analysis(self, dcc_args: DccExtrusionArgs, radii):
        sim_stats_pd = DccExtrusionRunner.read_stats_file(dcc_args.stats_file)
        for index, sim in sim_stats_pd.iterrows():
            dcc_args.output_folder = sim['sim_out_folder']
            if not os.path.exists(dcc_args.output_folder):
                logger.error(f'Missing {dcc_args.output_folder} so skip it!')
                continue
            self.sim_contact_radius_analysis(dcc_args, radii)

    def sim_contact_radius_analysis(self, dcc_args: DccExtrusionArgs, radii, dep_jobid=None):
        analysis_folders = []
        last_jobid = ''  # we'll need the last job to do the aggregation decay plot
        for rcp in radii:
            dcc_args_r = copy.deepcopy(dcc_args)
            cp, r = self.split_radius_contact_probability(rcp)
            dcc_args_r.radius_contact = r
            dcc_args_r.contact_probability = cp
            analysis_folders.append(dcc_args_r.default_analysis_folder())
            dcc_args_r.analyse = dcc_args_r.default_analysis_folder()
            last_jobid = self.run(dcc_args_r, stats_only=False, radii=radii, dep_jobid=dep_jobid)

            # do multi-decay plot on the done radius analysis
            last_jobid = self.run_multi_decay_plot(dcc_args, last_jobid)

        return last_jobid

    @staticmethod
    def split_radius_contact_probability(rcp):
        if rcp[-1].isdigit():
            r = float(rcp)
            cp = False
        else:
            r = float(rcp[:-1])
            cp = True
        return cp, r

    def run_multi_decay_plot(self, dcc_args: DccExtrusionArgs, dep_jobid):
        s_cmp_chrs = f'--cmp_chrs {" ".join(dcc_args.cmp_chrs)}' if dcc_args.cmp_chrs is not None else ''
        cmd_prefix = self.get_cmd_prefix(dcc_args, mpirun=False)
        cmd = f"{cmd_prefix} 3dpolys_le_runner multi_decay_plot -o {dcc_args.output_folder} -i {dcc_args.input_cfg} -e {dcc_args.exp_cool}" \
              f" {s_cmp_chrs}"
        return self._job_runner.run_cmd(cmd, dep_jobid)

    @staticmethod
    def multi_decay_plot(dcc_args: DccExtrusionArgs, replace=False):
        sub_folders = [dcc_args.analyse] if dcc_args.analyse else sorted([f.path for f in os.scandir(dcc_args.output_folder) if f.is_dir()])
        plots_folder = os.path.join(dcc_args.output_folder, ha.PLOTS_FOLDER)
        if plots_folder in sub_folders:
            sub_folders.remove(plots_folder)
        hic_h5_list = []
        for analysis_folder in sub_folders:  # list(sub_folders[1:2]):
            # find the last hic
            hic_h5_list.append(ha.get_last_hic(analysis_folder))
        hic_multi_decay_plot_log = ha.plot_distance_contact_prob_decay(hic_h5_list, exp_cool=dcc_args.exp_cool,
                                                                       output_folder=dcc_args.output_folder,
                                                                       res=dcc_args.resolution, confidence=0., replace=replace,
                                                                       chi2_mode=ha.CHI2_MODE_LOG)
        hic_multi_decay_plot_lin = ha.plot_distance_contact_prob_decay(hic_h5_list, exp_cool=dcc_args.exp_cool,
                                                                       output_folder=dcc_args.output_folder,
                                                                       res=dcc_args.resolution, confidence=0., replace=replace,
                                                                       chi2_mode=ha.CHI2_MODE_LINEAR)

    @staticmethod
    def multi_decay_exps_plot(dcc_args: DccExtrusionArgs, exp_cools, hic_chrs=None,
                              plot_format=ha.PLOT_FORMAT, replace=False):
        output_folder = dcc_args.output_folder
        if output_folder and not os.path.exists(output_folder):
            os.mkdir(output_folder)
        hic_multi_decay_plot_log = ha.plot_distance_contact_prob_decay(exp_cools[1:], hic_chrs=hic_chrs,
                                                                       exp_cool=exp_cools[0],
                                                                       tads=dcc_args.tads_boundary,
                                                                       output_folder=output_folder,
                                                                       res=dcc_args.resolution,
                                                                       confidence=0., replace=replace,
                                                                       chi2_mode=ha.CHI2_MODE_LOG,
                                                                       format=plot_format)
        hic_multi_decay_plot_lin = ha.plot_distance_contact_prob_decay(exp_cools[1:], hic_chrs=hic_chrs,
                                                                       exp_cool=exp_cools[0],
                                                                       tads=dcc_args.tads_boundary,
                                                                       output_folder=output_folder,
                                                                       res=dcc_args.resolution,
                                                                       confidence=0., replace=replace,
                                                                       chi2_mode=ha.CHI2_MODE_LINEAR,
                                                                       format=plot_format)

    @staticmethod
    def decay_plots(dcc_args: DccExtrusionArgs, res=ha.RESOLUTION, plot_format=ha.PLOT_FORMAT, replace=False,
                    use_threading=False):
        sim_stats_pd = DccExtrusionRunner.read_stats_file(dcc_args.stats_file)
        arg_exp_cool = dcc_args.exp_cool
        for index, sim in sim_stats_pd.iterrows():
            logger.info(
                f"{index}, hic:{sim['sim_hic_file']}, nlef:{sim['nlef']}, km:{sim['km']}, r:{sim['radius_contact']}")
            exp_cool = sim['exp_cool']
            hic_h5 = sim['sim_hic_file']
            if not os.path.exists(hic_h5):
                logger.error(f'Missing {hic_h5} so skip it!')
                continue

            # if contact_prob:
            if (arg_exp_cool and arg_exp_cool != EXP_COOL_AS_STATS) or not arg_exp_cool:  # also with no exp_cool
                exp_cool = arg_exp_cool
            else:
                exp_cool = sim['exp_cool']
            output_folder = sim['output_folder']

            if use_threading:
                t = threading.Thread(target=ha.plot_distance_contact_prob_decay, args=([hic_h5]),
                                     kwargs={'exp_cool': exp_cool,'output_folder': output_folder, 'res': res,
                                             'replace': replace,
                                             'chi2_mode': ha.CHI2_MODE_LOG,
                                             'format': plot_format})
                t.start()
            else:
                hic_decay_prob_plot = ha.plot_distance_contact_prob_decay([hic_h5], exp_cool=exp_cool,
                                                                          output_folder=output_folder,
                                                                          res=res,
                                                                          replace=replace,
                                                                          chi2_mode=ha.CHI2_MODE_LOG,
                                                                          format=plot_format)

    @staticmethod
    def chip_seq_plots(dcc_args: DccExtrusionArgs, correlation=ha.DEFAULT_CHIP_CORRELATION, replace=True):
        sim_stats_pd = DccExtrusionRunner.read_stats_file(dcc_args.stats_file)
        exp_chip = dcc_args.exp_chip
        for index, sim in sim_stats_pd.iterrows():
            hic_h5 = sim['sim_hic_file']
            if not os.path.exists(hic_h5):
                logger.error(f'Missing {hic_h5} so skip it!')
                continue
            dcc_args.output_folder = os.path.dirname(hic_h5)  # output_folder = sim['output_folder']  # maybe this way
            chip_out_file = os.path.join(dcc_args.output_folder, ha.CHIP_OUT)
            boundary = dcc_args.boundary
            resolution =dcc_args.resolution

            if not os.path.exists(chip_out_file):
                logger.error(f'Output simulation data file {chip_out_file} is missing so skip it!')
                continue

            ha.plot_chip_seq(chip_out_file, exp_chip, boundary, resolution, correlation,
                             plot=True, replace=replace)

    def grid_nlef_km(self, dcc_args: DccExtrusionArgs, nlef_list, km_list, radii, replace):
        for nlef in nlef_list:
            for km in km_list:     # for km in np.arange(2.7e-3 / 5, 2 * 2.7e-3, km_step):
                dcc_args_grid = copy.deepcopy(dcc_args)
                dcc_args_grid.nlef = nlef
                dcc_args_grid.km = km
                cur_jobid = self.run(dcc_args_grid, radii=radii, replace=replace)


def main():
    args = cli_parser().parse_args(sys.argv[1:])

    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)

    if args.cmp_chrs is not None:
        # TODO check if this is needed anymore
        logger.info(f'Overwriting the default hic_analysis.CHR_SYNONYMS: { ",".join(args.cmp_chrs)} '
                       f'with which HiCs will be compared!')
        ha.CHR_SYNONYMS = args.cmp_chrs

    job_runner = CfgJobRunner(input_cfg=args.input_cfg, cmd_run_file=args.cmd_run_file)

    dcc_args = DccExtrusionArgs(boundary=args.boundary, lef_loading_sites=args.lef_loading_sites,
                                input_cfg=args.input_cfg, tads_boundary=args.tads_boundary,
                                stats_file=args.stats_file, exp_cool=args.exp_cool, nlef=args.nlef, km=args.km,
                                radius_contact=args.radius_contact, contact_probability=args.contact_probability,
                                boundary_direction=args.boundary_direction,
                                z_loop=args.z_loop, unidirectional=args.unidirectional, init_mode=args.init_mode,
                                stats=args.stats, all_stats=args.all_stats,  # simultaneously=args.simultaneously,
                                cmp_chrs=args.cmp_chrs, resolution=args.resolution,
                                output_folder=args.output_folder, analyse=args.analysis_folder)
    dcc_run = DccExtrusionRunner(job_runner=job_runner)

    # plot_format = args.plot_format # TODO add args.plot_format
    plot_format = job_runner.get_property(profile='stats', name='plot_format')
    if not plot_format:
        plot_format = ha.PLOT_FORMAT  # default in case of ""

    if args.run_command == 'grid_nlef_km':
        dcc_run.grid_nlef_km(dcc_args=dcc_args, nlef_list=args.nlef_list, km_list=args.km_list, radii=args.list_contact_radii,
                             replace=args.replace)
    elif args.run_command == 'new_stats':
        dcc_run.analysis_stats(dcc_args, new_stats=True, exp_cool=args.exp_cool)
    elif args.run_command == 'decay_plots':
        dcc_run.decay_plots(dcc_args, res=args.resolution, plot_format=plot_format, replace=args.replace, use_threading=args.threading)
    elif args.run_command == 'chip_seq_plots':
        dcc_run.chip_seq_plots(dcc_args, correlation=args.correlation, replace=args.replace)
    elif args.run_command == 'contact_radius_analysis':
        dcc_run.contact_radius_analysis(dcc_args, radii=args.list_contact_radii)
    elif args.run_command == 'multi_decay_plot':
        dcc_run.multi_decay_plot(dcc_args, replace=args.replace)
    elif args.run_command == 'multi_decay_exps_plot':
        dcc_run.multi_decay_exps_plot(dcc_args, args.exp_cools, hic_chrs=args.hic_chrs,
                                      plot_format=plot_format, replace=args.replace)
    elif args.run_command == 'run':
        dcc_run.run(dcc_args, radii=args.list_contact_radii, replace=args.replace)


if __name__ == '__main__':
    main()
