#!/usr/bin/env python3
"""Compatibility wrapper for the packaged hgrid boundary implementation."""

from pathlib import Path
import sys

try:
    from .hgrid_preproc.boundaries import *  # noqa: F401,F403
    from .hgrid_preproc.boundaries import main
except ImportError:  # Support direct execution from a source checkout.
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from stofs3d_setup.ops.Grid.hgrid_preproc.boundaries import *  # noqa: F401,F403
    from stofs3d_setup.ops.Grid.hgrid_preproc.boundaries import main


if __name__ == "__main__":
    raise SystemExit(main())
