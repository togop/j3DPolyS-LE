#! /usr/bin/env python

import logging
import re
import subprocess
import configparser
import os.path

from abc import ABC, abstractmethod

# Initialization
logger = logging.getLogger(__name__)


class JobRunner(ABC):

    def run_cmd(self, cmd, dep_jobid, profile='') -> str:
        start_cmd = self._get_start_cmd(dep_jobid, profile)
        full_cmd = f"{start_cmd} {cmd}"

        jobid = ''
        if self._cmd_in_shell(profile):
            status, jobout = subprocess.getstatusoutput(full_cmd)
            logger.info(f" call: {full_cmd}\n\t {jobout}")
            if status == 0:
                jobid = self._get_jobid(jobout, profile)
                logger.info(f"JobID is: {jobid}")
            else:
                logger.error(f"Error submitting Job: {full_cmd}")
        elif self._cmd_in_stdout(profile):
            print(full_cmd)
        else:
            cmd_in_file = self._cmd_in_file(profile)
            if cmd_in_file:
                if not os.path.isfile(cmd_in_file):
                    with open(cmd_in_file, 'w') as f:
                        print('', file=f)
                with open(cmd_in_file, 'a') as f:
                    print(full_cmd, file=f)
        return jobid

    @abstractmethod
    def _get_start_cmd(self, dep_jobid, profile) -> str:
        pass

    @abstractmethod
    def _get_jobid(self, jobout, profile) -> str:
        pass

    @abstractmethod
    def _cmd_in_shell(self, profile) -> bool:
        pass

    @abstractmethod
    def _cmd_in_stdout(self, profile) -> bool:
        pass

    @abstractmethod
    def _cmd_in_file(self, profile) -> str:
        pass


class CfgJobRunner(JobRunner):
    input_cfg: str
    cmd_in_file: str
    _config = configparser.ConfigParser()

    def __init__(self, input_cfg, cmd_in_file):
        self.input_cfg = input_cfg
        self.cmd_in_file = cmd_in_file
        self._config.read(self.input_cfg)

    def _get_start_cmd(self, dep_jobid, profile) -> str:
        cmd_job_dependency = self._get_property(profile, 'cmd_job_dependency')
        cmd_dep = cmd_job_dependency.replace('{jobid}', dep_jobid) if dep_jobid else ''
        cmd_prefix = self._get_property(profile, 'cmd_prefix')
        return cmd_prefix.replace('{cmd_job_dependency}', cmd_dep)

    def _get_jobid(self, jobout, profile) -> str:
        jobid_re = self._get_property(profile, 'jobid_re')
        return re.search(jobid_re, jobout)[0]

    def _get_property(self, profile, name):
        value = ''
        if profile:
            try:
                value = self._config.get('job_runner_' + profile, name)
            except (configparser.NoOptionError, configparser.NoSectionError) as e:
                value = ''
        if not value:
            value = self._config.get('job_runner', name)
        return value

    def _cmd_in_shell(self, profile) -> bool:
        cmd_in = self._get_property(profile, 'cmd_in')
        return cmd_in == 'shell'

    def _cmd_in_stdout(self, profile) -> bool:
        cmd_in = self._get_property(profile, 'cmd_in')
        return cmd_in == 'stdout'

    def _cmd_in_file(self, profile) -> str:
        cmd_in = self._get_property(profile, 'cmd_in')
        if cmd_in.startswith('file:'):
            cmd_in_file = cmd_in.split(':')[1]
            return cmd_in_file.replace('{cmd_in_file}', self.cmd_in_file)
        return ''
