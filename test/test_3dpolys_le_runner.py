#! /usr/bin/env python

import logging
import subprocess
import os.path
import shutil
import time
import filecmp
import pytest
from py3dpolys_le import plot_hic

# Initialization
logger = logging.getLogger(__name__)


def is_file_created(file_path, timeout=300):  # max 15min
    wait = 0
    step = 1  # sec
    while not os.path.exists(file_path) and wait < timeout:
        time.sleep(step)
        wait += step
    return os.path.isfile(file_path)


def cmp_file_txt(expected_shel_script, shel_script):
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
        return False
    return True


def test_tads_shell_input():
    shel_script = "./run_test_tads_shell_input.sh"
    if os.path.exists(shel_script):
        os.remove(shel_script)
    if os.path.isdir('out/test3'):
        shutil.rmtree('out/test3')
    plot_format = 'svg'
    cmd = f"3dpolys_le_runner run -o out/test3 --km 0.0017 --nlef 100 -bd 1 -z -u -im h -lr 3.55 7.1 " \
          f"-i ./test/test_tads_shell_input.cfg --cmd_run_file {shel_script}"
    print(f"call: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    assert is_file_created(shel_script), f"Shell script {shel_script} file is not created"
    expected_shel_script = "./test/expected/run_test_tads_shell_input.sh"
#TODO    assert cmp_file_txt(expected_shel_script, shel_script), f"Shell script {shel_script} is not as expected: "
    subprocess.run(f"chmod +x {shel_script}", shell=True, check=True)
    subprocess.run(shel_script, shell=True, check=True)
    assert is_file_created(f'out/test3/r3.55/hic_001_cool.{plot_format}'),\
        f"File out/test3/r3.55/hic_001_cool.{plot_format} is missing"
    assert is_file_created(f'out/test3/r3.55/hic_002_cool.{plot_format}'),\
        f"File out/test3/r3.55/hic_002_cool.{plot_format} is missing"
    assert is_file_created(f'out/test3/r7.10/hic_001_cool.{plot_format}'), \
        f"File out/test3/r7.10/hic_001_cool.{plot_format} is missing"
    assert is_file_created(f'out/test3/r7.10/hic_002_cool.{plot_format}'), \
        f"File out/test3/r7.10/hic_002_cool.{plot_format} is missing"


def test_no_tads_shell():
    shel_script = "./run_test_no_tads_shell_input.sh"
    if os.path.exists(shel_script):
        os.remove(shel_script)
    cmd = f"3dpolys_le_runner run -i ./test/test_no_tads_shell_input.cfg " \
          f"--cmd_run_file {shel_script}"
    print(f"call: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    assert is_file_created(shel_script), f"Shell script {shel_script} file is not created"
    # assert cmp_file_txt(expected_shel_script, shel_script), f"Shell script {shel_script} is not as expected: "
    subprocess.run(f"chmod +x {shel_script}", shell=True, check=True)
    subprocess.run(shel_script, shell=True, check=True)


def test_shell_container():
    shel_script = "./run_test_shell_container_input.sh"
    if os.path.exists(shel_script):
        os.remove(shel_script)
    if os.path.isdir('out/test_shell_container'):
        shutil.rmtree('out/test_shell_container')
    os.mkdir('out/test_shell_container')  # create new
    cmd = "3dpolys_le_runner run -i ./test/test_shell_container_input.cfg -o out/test_shell_container " \
          f"--cmd_run_file {shel_script}"
    print(f"call: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    assert is_file_created(shel_script), f"Shell script {shel_script} file is not created"
    expected_shel_script = "./test/expected/run_test_shell_container_input.sh"
    assert cmp_file_txt(expected_shel_script, shel_script), f"Shell script {shel_script} is not as expected: "

def test_tads_init_continue():
    # init
    if os.path.isdir('out/test_init'):
        shutil.rmtree('out/test_init')
    cmd = f"mpirun 3dpolys_le -o:out/test_init ./test/test_tads_init.cfg"
    print(f"call: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    assert is_file_created(f'out/test_init/3dpoys_le.cfg'), \
        f"File out/test_init/3dpoys_le.cfg is missing"
    assert is_file_created(f'out/test_init/config.out'), \
        f"File out/test_init/config.out is missing"
    assert is_file_created(f'out/test_init/contact.out'), \
        f"File out/test_init/contact.out is missing"
    # continue
    if os.path.isdir('out/test_continue'):
        shutil.rmtree('out/test_continue')
    cmd = f"mpirun 3dpolys_le -o:out/test_continue -im:s=out/test_init ./test/test_tads_continue.cfg"
    print(f"call: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    assert is_file_created(f'out/test_continue/3dpoys_le.cfg'), \
        f"File out/test_continue/3dpoys_le.cfg is missing"
    assert is_file_created(f'out/test_continue/config.out'), \
        f"File out/test_continue/config.out is missing"
    assert is_file_created(f'out/test_continue/contact.out'), \
        f"File out/test_continue/contact.out is missing"


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)

    test_tads_shell_input()
    test_no_tads_shell()
    test_shell_container()
    test_tads_init_continue()

