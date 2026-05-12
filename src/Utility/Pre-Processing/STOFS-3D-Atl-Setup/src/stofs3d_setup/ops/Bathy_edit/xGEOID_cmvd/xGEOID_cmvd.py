import numpy as np
from copy import deepcopy
from coastalmodeling_vdatum import vdatum


def convert2xgeoid_cmvd(hgrid_obj_lonlat_navd88):
    """
    Convert the hgrid depth (positive down) from NAVD88 to xGEOID using Felicio's alternative tool cmvd:
    https://github.com/oceanmodeling/coastalmodeling-vdatum

    hgrid_obj_lonlat_navd88: Hgrid object, the depth is expected to be in NAVD88, and will be converted to xGEOID.
            hgrid_obj_lonlat_navd88.x and hgrid_obj_lonlat_navd88.y are expected to be in lat/lon
    """
    hg = deepcopy(hgrid_obj_lonlat_navd88)

    # coastalmodeling-vdatum expects positive z overland and negative z under water, thus multiple by -1
    z_in = -hgrid_obj_lonlat_navd88.dp
    x, y, z_out = vdatum.convert("navd88", "xgeoid20b", hg.y, hg.x, z_in, online=True, epoch=None)
    z_out = -z_out  # convert back to positive down
    # if the conversion fails, vdatum returns inf, we keep the original z in this case
    valid = np.isfinite(z_out)
    hg.dp[valid] = z_out[valid]

    return hg
