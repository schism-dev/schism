#!/usr/bin/env python3

import geopandas as gpd
import numpy as np
from pathlib import Path

from matplotlib import pyplot as plt

from pylib_essentials.schism_file import read_schism_hgrid_cached, grd2sms, schism_grid


def load_NCF(hgrid_obj:schism_grid, NCF_shpfile:Path, expansion:float=4.0):
    '''
    Load the maintained depth from National Channel Framework data into the hgrid.
    The original NCF polygons are enlarged to accommodate the mismatch between the hgrid and the NCF data.
    The maintained depth is converted from feet to meters.
    '''

    print('Loading NCF data...\n')
    NCF_data = gpd.read_file(NCF_shpfile)
    # remove the polygons that are not in the bounding box of the hgrid
    NCF_data = NCF_data.cx[hgrid_obj.x.min():hgrid_obj.x.max(), hgrid_obj.y.min():hgrid_obj.y.max()]
    # convert any multipolygons to polygons
    NCF_data = NCF_data.explode(index_parts=True)

    print('Enlarging the NCF polygons...\n')
    NCF_data['geometry'] = NCF_data['geometry'].to_crs('esri:102008').buffer(expansion).to_crs('epsg:4326')

    print('setting dp at points inside the NCF polygons...\n')
    # put hgrid points into a Point GeoDataFrame
    hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(hgrid_obj.x, hgrid_obj.y), crs='epsg:4326')
    # determine which points are inside the polygons and get the maintained depth
    joined_gdf = gpd.sjoin(hg_points, NCF_data, how="inner", predicate='within')
    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons
    dp_NCF = np.ones_like(hgrid_obj.dp, dtype=float) * -9999  # initialize with a large negative number
    dp_NCF[idx] = joined_gdf['depthmaint'].to_numpy() * 0.3048  # convert from feet to meters

    # diagnostic plot
    # plt.figure()
    # hgrid_obj.plot(value=dp_NCF.astype(int), fmt=1)
    # plt.show()

    hgrid_obj.dp = dp_NCF

    return hgrid_obj

def plot_diagnostic(hgrid_obj:schism_grid, dp:np.ndarray):
    # # plot to check the changed nodes
    # plt.figure()
    # idx = abs(dp_orig - hgrid_obj.dp) > 1e-5
    # # add a 1:1 line for reference from the lower left to upper right
    # # Determine the range for the 1:1 line
    # min_val = min(min(dp_orig[idx]), min(hgrid_obj.dp[idx]))
    # max_val = max(max(dp_orig[idx]), max(hgrid_obj.dp[idx]))
    # # Add a 1:1 line for reference
    # plt.plot([min_val, max_val], [min_val, max_val], 'k--')
    # plt.scatter(dp_orig[idx], hgrid_obj.dp[idx])
    # plt.axis('equal')
    # plt.show()
    pass

if __name__ == '__main__':
    # sample usage
    # -----------------inputs-----------------------------
    hg_file = Path('/sciclone/schism10/Hgrid_projects/TMP/DEM_edit/NCF/hgrid_dem_levee_loaded.gr3')
    # -----------------end inputs-----------------------------

    hgrid_obj = schism_grid(str(hg_file))
    dp_orig = hgrid_obj.dp.copy()

    # the NCF data already defines the region
    hgrid_obj = load_NCF(hgrid_obj=hgrid_obj, NCF_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/channel_quarter_NCF.shp'))
    hgrid_obj.dp = np.maximum(dp_orig, hgrid_obj.dp)  # change the depth only if it is deeper than the original depth
    grd2sms(hgrid_obj, (f"{hg_file.parent}/{hg_file.stem}_NCF_loaded.2dm"))
    hgrid_obj.save(f"{hg_file.parent}/{hg_file.stem}_NCF_loaded.gr3", fmt=1)

    pass
