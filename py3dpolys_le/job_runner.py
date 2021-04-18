#! /usr/bin/env python

import logging
import re
import subprocess

from abc import ABC, abstractmethod

# Initialization
logger = logging.getLogger(__name__)


class JobRunner(ABC):

    def run_cmd(self, cmd, dep_jobid) -> str:
        start_cmd = self._get_start_cmd(dep_jobid)
        slurm_cmd = f"{start_cmd} {cmd}"

        jobid = ''
        status, jobout = subprocess.getstatusoutput(slurm_cmd)
        logger.info(f" call: {slurm_cmd}\n\t {jobout}")
        if status == 0:
            jobid = self._get_jobid(jobout)
            logger.info(f"JobID is: {jobid}")
        else:
            logger.error(f"Error submitting Job: {slurm_cmd}")
        return jobid

    @abstractmethod
    def _get_start_cmd(self, dep_jobid) -> str:
        pass

    @abstractmethod
    def _get_jobid(self, jobout) -> str:
        pass


class ShellJobRunner(JobRunner):
    def _get_start_cmd(self, dep_jobid) -> str:
        return ''

    def _get_jobid(self, jobout) -> str:
        return ''


class SlurmJobRunner(JobRunner):
    def _get_start_cmd(self, dep_jobid) -> str:
        return f"sbatch --dependency=afterany:{dep_jobid}" if len(dep_jobid) != 0 and not dep_jobid.isspace() else "sbatch"

    def _get_jobid(self, jobout) -> str:
        return re.findall(r"\d+$", jobout)[0]


class DummyJobRunner(JobRunner):

    def _get_start_cmd(self, dep_jobid) -> str:
        return f"echo depends={dep_jobid} cmd:"

    def _get_jobid(self, jobout):
        return str(abs(hash(jobout)))
