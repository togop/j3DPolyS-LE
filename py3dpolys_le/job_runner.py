#! /usr/bin/env python

import logging
import re
import subprocess
import configparser

from abc import ABC, abstractmethod

# Initialization
logger = logging.getLogger(__name__)


class JobRunner(ABC):

    def run_cmd(self, cmd, dep_jobid, profile='') -> str:
        start_cmd = self._get_start_cmd(dep_jobid, profile)
        full_cmd = f"{start_cmd} {cmd}"

        jobid = ''
        status, jobout = subprocess.getstatusoutput(full_cmd)
        logger.info(f" call: {full_cmd}\n\t {jobout}")
        if status == 0:
            jobid = self._get_jobid(jobout)
            logger.info(f"JobID is: {jobid}")
        else:
            logger.error(f"Error submitting Job: {full_cmd}")
        return jobid

    @abstractmethod
    def _get_start_cmd(self, dep_jobid, profile) -> str:
        pass

    @abstractmethod
    def _get_jobid(self, jobout) -> str:
        pass


class CfgJobRunner(JobRunner):
    cmd_prefix: str
    cmd_prefix_sim: str
    cmd_prefix_analysis: str
    cmd_prefix_stats: str
    jobid_re: str
    cmd_job_dependency: str

    def __init__(self, input_cfg):
        config = configparser.ConfigParser()
        config.read(input_cfg)
        self.cmd_prefix = config.get('job_runner', 'cmd_prefix')
        self.cmd_prefix_sim = config.get('job_runner', 'cmd_prefix_sim')
        self.cmd_prefix_analysis = config.get('job_runner', 'cmd_prefix_analysis')
        self.cmd_prefix_stats = config.get('job_runner', 'cmd_prefix_stats')
        self.jobid_re = config.get('job_runner', 'jobid_re')
        self.cmd_job_dependency = config.get('job_runner', 'cmd_job_dependency')

    def _get_start_cmd(self, dep_jobid, profile) -> str:
        if profile == 'sim':
            cmd_dep = self.cmd_prefix_sim.replace('{jobid}', dep_jobid)
        elif profile == 'analysis':
            cmd_dep = self.cmd_prefix_analysis.replace('{jobid}', dep_jobid)
        elif profile == 'stats':
            cmd_dep = self.cmd_prefix_stats.replace('{jobid}', dep_jobid)
        else:  # various plots
            cmd_dep = self.cmd_job_dependency.replace('{jobid}', dep_jobid)
        return self.cmd_prefix.replace('{cmd_job_dependency}', cmd_dep)

    def _get_jobid(self, jobout) -> str:
        return re.search(self.jobid_re, jobout)[0]
