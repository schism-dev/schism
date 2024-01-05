import geopandas as gpd
import numpy as np
from pathlib import Path

from matplotlib import pyplot as plt

from pylib_essentials.schism_file import read_schism_hgrid_cached, grd2sms, schism_grid


def load_NCF(hg:schism_grid, NCF_shpfile:Path):
    '''
    Load the maintained depth from National Channel Framework data into the hgrid.
    The original NCF polygons are enlarged to accommodate the mismatch between the hgrid and the NCF data.
    The maintained depth is converted from feet to meters.
    '''

    # get the nodes in the NCF (National Channel Framework) polygons
    NCF_shpfile = Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/channel_quarter_NCF.shp')
    NCF_data = gpd.read_file(NCF_shpfile)

    # remove the polygons that are not in the bounding box of the hgrid
    NCF_data = NCF_data.cx[hg.x.min():hg.x.max(), hg.y.min():hg.y.max()]
    # convert any multipolygons to polygons
    NCF_data = NCF_data.explode()

    # enlarge the polygons by 4 m
    NCF_data['geometry'] = NCF_data['geometry'].to_crs('esri:102008').buffer(4).to_crs('epsg:4326')
    # extract the polygons to a list of nx2 numpy arrays
    # NCF_polygons = [np.array(poly.exterior.coords) for poly in NCF_data['geometry']]

    # put hgrid points into a Point GeoDataFrame
    hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(hg.x, hg.y), crs='epsg:4326')
    # determine which points are inside the polygons and get the maintained depth
    joined_gdf = gpd.sjoin(hg_points, NCF_data, how="inner", predicate='within')
    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons
    dp_NCF = np.ones_like(hg.dp, dtype=float) * -9999  # initialize with a large negative number
    dp_NCF[idx] = joined_gdf['depthmaint'].to_numpy() * 0.3048  # convert from feet to meters

    # diagnostic plot
    # plt.figure()
    # hg.plot(value=dp_NCF.astype(int), fmt=1)
    # plt.show()

    return dp_NCF

def plot_diagnostic(hg:schism_grid, dp:np.ndarray):
    # # plot to check the changed nodes
    # plt.figure()
    # idx = abs(dp_orig - hg.dp) > 1e-5
    # # add a 1:1 line for reference from the lower left to upper right
    # # Determine the range for the 1:1 line
    # min_val = min(min(dp_orig[idx]), min(hg.dp[idx]))
    # max_val = max(max(dp_orig[idx]), max(hg.dp[idx]))
    # # Add a 1:1 line for reference
    # plt.plot([min_val, max_val], [min_val, max_val], 'k--')
    # plt.scatter(dp_orig[idx], hg.dp[idx])
    # plt.axis('equal')
    # plt.show()
    pass

if __name__ == '__main__':
    # -----------------inputs-----------------------------
    hg_file = Path('./hgrid_dem_levee_loaded.gr3')
    # -----------------end inputs-----------------------------

    hg = schism_grid(str(hg_file))
    dp_orig = hg.dp.copy()

    # the NCF data already defines the region
    dp_from_NCF = load_NCF(hg=hg, NCF_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/channel_quarter_NCF.shp'))
    hg.dp = np.maximum(dp_orig, dp_from_NCF)  # change the depth only if it is deeper than the original depth
    grd2sms(hg, (f"{hg_file.parent}/{hg_file.stem}_NCF_loaded.2dm"))
    hg.save(f"{hg_file.parent}/{hg_file.stem}_NCF_loaded.gr3", fmt=1)

    pass
