import numpy as np
from scipy.interpolate import griddata

from pylib import read, read_schism_reg, inside_polygon


def interpolate_upper_hudson(
    hgrid_obj,
    reference_hgrid_file,
    region_file,
):
    """
    Interpolate bathymetry from a reference hgrid onto the base hgrid
    within the specified region.
    """

    ref = read(str(reference_hgrid_file))
    bp = read_schism_reg(str(region_file))

    # base-grid nodes inside Upper Hudson region
    sind = inside_polygon(
        np.c_[hgrid_obj.x, hgrid_obj.y],
        bp.x,
        bp.y,
    ) == 1

    print(f'Upper Hudson interpolation: {sind.sum()} nodes')

    # interpolate reference bathymetry to base-grid node locations
    dp_interp = griddata(
        np.c_[ref.x, ref.y],
        ref.dp,
        np.c_[hgrid_obj.x[sind], hgrid_obj.y[sind]],
        method='linear',
    )

    # fallback for nodes outside interpolation convex hull
    nan_idx = np.isnan(dp_interp)

    if np.any(nan_idx):
        dp_interp[nan_idx] = griddata(
            np.c_[ref.x, ref.y],
            ref.dp,
            np.c_[
                hgrid_obj.x[sind][nan_idx],
                hgrid_obj.y[sind][nan_idx],
            ],
            method='nearest',
        )

    hgrid_obj.dp[sind] = dp_interp

    return hgrid_obj
