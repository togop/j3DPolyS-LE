#! /usr/bin/env python

import logging
import subprocess
import os.path
import time
import filecmp
import pytest
from py3dpolys_le import plot_hic

from abc import ABC, abstractmethod

# Initialization
logger = logging.getLogger(__name__)


def is_file_created(file_path, timeout=300):  # max 15min
    wait = 0
    step = 1  # sec
    while not os.path.exists(file_path) and wait < timeout:
        time.sleep(step)
        wait += step
    return os.path.isfile(file_path)


def test_no_tads_shell():
    cmd = "3dpolys_le_runner run -i ./test/test_no_tads_shell_input.cfg " \
          "--cmd_run_file ./run_test_no_tads_shell_input.sh"
    print(f"call: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    shel_script = "./run_test_no_tads_shell_input.sh"
    assert is_file_created(shel_script), f"Shell script {shel_script} file is not created"
    subprocess.run(f"chmod +x {shel_script}", shell=True, check=True)
    subprocess.run(shel_script, shell=True, check=True)


def test_shell_container():
    shel_script = "./run_test_shell_container_input.sh"
    if os.path.exists(shel_script):
        os.remove(shel_script)
    cmd = "3dpolys_le_runner run -i ./test/test_shell_container_input.cfg -o out_test_shell_container " \
          f"--cmd_run_file {shel_script}"
    print(f"call: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    assert is_file_created(shel_script), f"Shell script {shel_script} file is not created"
    expected_shel_script = "./test/expected/run_test_shell_container_input.sh"
    status, stout = subprocess.getstatusoutput(f'./test/test_cmp.sh {expected_shel_script} {shel_script}')
    if status != 0:  # stout.strip().index('TEST OK:') < 0:  # not filecmp.cmp(shel_script, expected_shel_script):
        f = open(expected_shel_script, 'r')
        content = f.read()
        print('Expected:\n')
        print(content)
        f.close()

        f = open(shel_script, 'r')
        content = f.read()
        print('Got:\n')
        print(content)
        f.close()
        assert False, f"Shell script {shel_script} is not as expected: "


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)

    test_shell_container()
    #test_no_tads_shell()

