import os
import pkg_resources
# import subprocess
import sys


def main():
    bin_3dpolys_le = pkg_resources.resource_filename(__name__, 'bin/3dpolys_le')
    os.system(bin_3dpolys_le + " " + " ".join(sys.argv[1:]))
    # status, jobout = subprocess.getstatusoutput(bin_3dpolys_le + " ".join(sys.argv[1:]))
    # print(jobout)


if __name__ == '__main__':
    main()
