"""
This script provides utility methods used by RiverMapper
"""


import os
import shutil


def silentremove(filenames):
    if isinstance(filenames, str):
        filenames = [filenames]

    for filename in filenames:
        try:
            os.remove(filename)
        except IsADirectoryError:
            shutil.rmtree(filename)
        except FileNotFoundError:
            pass  # missing_ok
