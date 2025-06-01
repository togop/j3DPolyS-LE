#! /usr/bin/env python

import logging
import re
import subprocess
import configparser
import os.path

from abc import ABC, abstractmethod

# Initialization
logger = logging.getLogger(__name__)

CFG_SECTION_3DPOLYS_LE = '3dpolys_le'


class JobRunner(ABC):

    def run_cmd(self, cmd, dep_jobid, profile='') -> str:
        start_cmd = self._get_start_cmd(dep_jobid, profile)
        full_cmd = f"{start_cmd} {cmd}"

        jobid = ''
        if self._cmd_run_shell(profile):
            status, jobout = subprocess.getstatusoutput(full_cmd)
            logger.info(f" call: {full_cmd}\n\t {jobout}")
            if status == 0:
                jobid = self._get_jobid(jobout, profile)
                logger.info(f"JobID is: {jobid}")
            else:
                logger.error(f"Error submitting Job: {full_cmd}")
        elif self._cmd_run_stdout(profile):
            print(full_cmd)
        else:
            cmd_run_file = self._cmd_run_file(profile)
            if cmd_run_file:
                if not os.path.isfile(cmd_run_file):
                    with open(cmd_run_file, 'w') as f:
                        print('#! /bin/bash', file=f)  # TODO make ! /bin/bash directive parameterized
                with open(cmd_run_file, 'a') as f:
                    print(full_cmd, file=f)
        return jobid

    @abstractmethod
    def _get_start_cmd(self, dep_jobid, profile) -> str:
        pass

    @abstractmethod
    def _get_jobid(self, jobout, profile) -> str:
        pass

    @abstractmethod
    def _cmd_run_shell(self, profile) -> bool:
        pass

    @abstractmethod
    def _cmd_run_stdout(self, profile) -> bool:
        pass

    @abstractmethod
    def _cmd_run_file(self, profile) -> str:
        pass


class CfgJobRunner(JobRunner):
    input_cfg: str
    cmd_run_file: str
    _config = configparser.ConfigParser()

    def __init__(self, input_cfg, cmd_run_file=None):
        self.input_cfg = input_cfg
        self.cmd_run_file = cmd_run_file
        self._config.read(self.input_cfg)

    def _get_start_cmd(self, dep_jobid, profile) -> str:
        cmd_job_dependency = self.get_property(profile, 'cmd_job_dependency')
        cmd_dep = cmd_job_dependency.replace('{jobid}', dep_jobid) if dep_jobid else ''
        cmd_prefix = self.get_property(profile, 'cmd_prefix')
        return cmd_prefix.replace('{cmd_job_dependency}', cmd_dep)

    def _get_jobid(self, jobout, profile) -> str:
        jobid_re = self.get_property(profile, 'jobid_re')
        return re.search(jobid_re, jobout)[0]

    def get_sim_property(self, name):
        value = ''
        try:
            value = self._config.get(CFG_SECTION_3DPOLYS_LE, name)
        except (configparser.NoOptionError, configparser.NoSectionError) as e:
            value = ''
        return value

    def get_property(self, profile, name, default=''):
        value = default
        if profile:
            try:
                value = self._config.get('job_runner_' + profile, name)
            except (configparser.NoOptionError, configparser.NoSectionError) as e:
                value = default
        if not value:
            try:
                value = self._config.get('job_runner', name)
            except (configparser.NoOptionError, configparser.NoSectionError) as e:
                value = default
        return value

    def _cmd_run_shell(self, profile) -> bool:
        cmd_run = self.get_property(profile, 'cmd_run')
        return cmd_run == 'shell' and not self.cmd_run_file

    def _cmd_run_stdout(self, profile) -> bool:
        cmd_run = self.get_property(profile, 'cmd_run')
        return cmd_run == 'stdout'  and not self.cmd_run_file

    def _cmd_run_file(self, profile) -> str:
        if self.cmd_run_file:
            return self.cmd_run_file
        cmd_run = self.get_property(profile, 'cmd_run')
        if cmd_run.startswith('file:'):
            cmd_run_file = cmd_run.split(':')[1]
            return cmd_run_file  #  cmd_run_file.replace('{cmd_run_file}', self.cmd_run_file)
        return ''


def convert_dat_to_cfg(input_dat):
    _config = configparser.ConfigParser()
    _config.add_section(CFG_SECTION_3DPOLYS_LE)
    with open(input_dat) as fp:
        while True:
            line = fp.readline()
            if not line:
                break

            val, name = [x.strip() for x in line.rsplit("::")]
            _config.set(CFG_SECTION_3DPOLYS_LE, name, val)

    input_cfg = ".cfg".join(input_dat.rsplit(".dat", 1))
    with open(input_cfg, 'w') as cf:
        _config.write(cf)
    return input_cfg
