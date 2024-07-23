"""logging configuration for STOFS-3D-ATL."""

import logging
import os


def setup_logging(level=logging.DEBUG):
    '''Setup logging for STOFS-3D-ATL.'''
    logging.basicConfig(
        format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
        level=level,
        stream=os.sys.stdout
    )
