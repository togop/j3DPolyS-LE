import os
from importlib.resources import files, as_file
import sys


def main():
    with as_file(files('py3dpolys_le.bin') / '3dpolys_le') as bin_3dpolys_le_path:
        bin_3dpolys_le = str(bin_3dpolys_le_path)
    #    bin_3dpolys_le = resources.path('py3dpolys_le.bin', '3dpolys_le')
        cmd = str(bin_3dpolys_le) + " " + " ".join(sys.argv[1:])
        print(f'cmd:{cmd}')
        os.system(cmd)
    # status, jobout = subprocess.getstatusoutput(bin_3dpolys_le + " ".join(sys.argv[1:]))
    # print(jobout)


if __name__ == '__main__':
    main()
