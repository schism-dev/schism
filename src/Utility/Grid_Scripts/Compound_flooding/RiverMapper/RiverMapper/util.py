"""
This script provides utility methods used by RiverMapper
"""


import os
import errno
import shutil


def silentremove(filename):
    try:
        os.remove(filename)
    except IsADirectoryError:
        shutil.rmtree(filename)
    except FileNotFoundError:
        pass  # missing_ok
