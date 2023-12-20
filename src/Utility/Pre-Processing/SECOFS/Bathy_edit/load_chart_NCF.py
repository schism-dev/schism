import geopandas as gpd
import numpy as np
from pathlib import Path

from scipy.spatial import KDTree
from matplotlib import pyplot as plt

from pylib_essentials.schism_file import read_schism_hgrid_cached, grd2sms, schism_grid
from pylib_essentials.utility_functions import inside_polygon

def prep_chart_data(sounding_shpfile:Path):
    '''
    Prepare the chart data for loading into the hgrid.
    The original chart data stores coordinates in the geometry column.
    This function extracts the (x, y) to form a new geometry,
    and adds the z coordinate as a column so that it can be better visualized in GIS.

    Visualization is needed to check if the chart data is dense enough for interpolation.
    If not, the chart data needs to be manually edited to add more points.
    '''

    # sounding_shpfile = Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/savannah_cooper_sounding_2.shp')

    sounding_data = gpd.read_file(sounding_shpfile)
    # only keep the geometry column
    sounding_data = sounding_data['geometry']
    # extract xyz coordinates into an array
    sounding_xyz = np.array([point.coords[0] for point in sounding_data])

    # create a new geodataframe with the xyz coordinates, and set z as a column
    sounding_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(sounding_xyz[:,0], sounding_xyz[:,1]))
    sounding_gdf['z'] = sounding_xyz[:,2]

    # write gdf to shapefile
    output_file = f"{sounding_shpfile.parent}/{sounding_shpfile.stem}_xyz{sounding_shpfile.suffix}"
    sounding_gdf.to_file(output_file)

    return output_file

def load_chart(hg:schism_grid, sounding_shpfile:Path, region_shpfile:Path, crs_region:str='esri:102008'):
    '''
    Load the chart data into the hgrid.
    The sounding_shpfile is the chart data in xyz format, which is ususally based on MLLW.
    Convert the chart data to hgrid's datum if needed.
    A region file is needed to limit the hgrid nodes to be loaded.
    The regions are usually manually defined in SMS, which assumes esri:102008.
    The regions ususally need to be shrinked, especially when they are made from existing SMS channel polygons,
    so that the bank nodes are not selected and only the nodes in the channel are selected.
    '''

    # read the sounding data
    sounding_data = gpd.read_file(sounding_shpfile)
    # extract xy from geometry
    xy = np.array([point.coords[0] for point in sounding_data['geometry']])
    sounding_xyz = np.c_[xy, sounding_data['z'].values]

    # prepare the polygon data
    channel_shpfile = Path(region_shpfile, crs=crs_region)
    channel_polys = gpd.read_file(channel_shpfile)

    # shrink the polygons to exclude the bank nodes
    channel_polys['geometry'] = channel_polys['geometry'].buffer(-1)

    # intersect hgrid points with the polygons
    hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(hg.x, hg.y), crs='epsg:4326').to_crs(crs_region)
    joined_gdf = gpd.sjoin(hg_points, channel_polys, how="inner", predicate='within')

    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons
    in_channel = np.zeros_like(hg.dp, dtype=bool)
    in_channel[idx] = True

    # diagnostic plot
    plt.figure()
    hg.plot(value=in_channel.astype(int), fmt=1)

    # for inpolygon nodes, find the nearest sounding_xyz point using k-d tree
    idx = KDTree(sounding_xyz[:, :2]).query(np.c_[hg.x, hg.y])[1]
    dp_sounding = sounding_xyz[idx, 2]

    return dp_sounding


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
    joined_gdf = gpd.sjoin(hg_points, NCF_data, how="inner", op='within')
    idx = joined_gdf.index.to_numpy()  # get the indices of the points inside the polygons
    dp_NCF = np.ones_like(hg.dp, dtype=float) * -9999  # initialize with a large negative number
    dp_NCF[idx] = joined_gdf['depthmaint'].to_numpy() * 0.3048  # convert from feet to meters

    # diagnostic plot
    plt.figure()
    hg.plot(value=dp_NCF.astype(int), fmt=1)

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
    hg_file = Path('/sciclone/schism10/feiye/Test/RUN02b_JZ/hgrid.gr3')
    hg = read_schism_hgrid_cached(hg_file)

    # the chart data needs a manual region file to limit the nodes to be loaded
    dp_from_chart = load_chart(
        hg=hg,
        sounding_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/savannah_cooper_sounding_2_xyz_edited.shp'),
        region_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/manual_chart_region.shp'), crs_region='esri:102008'
    )
    hg.dp = np.maximum(hg.dp, dp_from_chart)  # change the depth only if it is deeper than the original depth

    # the NCF data already defines the region
    dp_from_NCF = load_NCF(hg=hg, NCF_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/channel_quarter_NCF.shp'))
    hg.dp = np.maximum(hg.dp, dp_from_NCF)  # change the depth only if it is deeper than the original depth

    grd2sms(hg, (f"{hg_file.parent}/{hg_file.stem}_chart_NCF_loaded.2dm"))

    pass