"""Paths to shared STOFS-3D-Atlantic v7.4 workflow input data."""

from pathlib import Path


WORKFLOW_CONSTANTS_DIR = Path(
    "/sciclone/schism10/Hgrid_projects/WORKFLOW_CONSTANTS_stofs3d_v7.4"
)
NWM_ECGC_SHAPEFILE = WORKFLOW_CONSTANTS_DIR / "NWM" / "ecgc" / "ecgc.shp"
