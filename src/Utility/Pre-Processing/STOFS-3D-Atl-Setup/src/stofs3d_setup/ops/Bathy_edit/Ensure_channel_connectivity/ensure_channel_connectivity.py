import os
import numpy as np
from copy import deepcopy
import geopandas as gpd
from RiverMapper.rivers import Rivers
from RiverMapper.SMS import SMS_MAP
from sklearn.neighbors import KDTree


def dredge_river_transects(
    rivers: Rivers,
    region_gdf: gpd.GeoDataFrame = None,
    hgrid_obj=None,  # schism_hgrid object read by the read() function in pylib
    min_channel_depth=1.0,
    measured_from_high_bank=True,
    output_dir='./'
):
    '''
    Dredge the inner arcs of each river transect to maintain a minimum elevation drop
    from bank to thalweg, thus maintaining channel connectivity

    Inputs:
    - rivers: Rivers object from reading the RiverMapper diagnostic output file
    - region_gdf: gpd.GeoDataFrame, region of interest, in which the dredging is performed.
        It must have a coordinate reference system (CRS) defined.
    - hgrid_obj: schism_hgrid object with bathymetry loaded, assuming lon/lat
    - min_channel_depth: float, minimum channel depth to dredge.
        The depth is measured from the higher bank elevation to an inner arc node.
    - measured_from_high_bank: bool, whether the channel depth is measured from the higher bank elevation.
    - output_dir: str, directory to save the dredged mesh and diagnostic files
    '''

    print('getting river arcs z from the mesh ...')
    rivers.mesh_dp2riverarc_z(hgrid_obj)

    print('dredging river transects ...')
    dredged_points = rivers.dredge_inner_arcs(
        region_gdf=region_gdf, min_channel_depth=min_channel_depth,
        inner_most_dredge=False,  # dredge all inner arcs, won't work if outer arcs are present
        measured_from_high_bank=measured_from_high_bank
    )

    print('mapping dredged points to the mesh ...')
    _, idx = KDTree(np.c_[hgrid_obj.x, hgrid_obj.y]).query(dredged_points[:, :2])
    # update the mesh
    hgrid_dredged = deepcopy(hgrid_obj)
    hgrid_dredged.dp[np.squeeze(idx)] = np.maximum(
        hgrid_dredged.dp[np.squeeze(idx)], dredged_points[:, 2]
    )

    print('saving dredged mesh ...')
    os.makedirs(output_dir, exist_ok=True)
    hgrid_dredged.grd2sms(output_dir + '/hgrid_dredged.2dm')  # SMS format
    hgrid_dredged.save(output_dir + '/hgrid_dredged.gr3', fmt=1)  # SCHISM format

    return hgrid_dredged


def ensure_channel_connectivity(
    hgrid_obj, min_channel_depth=1.0, measured_from_high_bank=True,
    river_extra_info_map_file=None,
    region_gdf_file=None, exclude_region_gdf_file_list=None,
    output_dir=None
):
    '''
    Ensure channel connectivity by dredging river transects defined by RiverMapper.
    The river transects are dredged to maintain a minimum elevation drop from bank to thalweg,
    thus maintaining channel connectivity.

    Inputs:
    - hgrid_obj: schism_hgrid object with bathymetry loaded, assuming lon/lat
    - river_extra_info_map_file: str, path to the RiverMapper diagnostic output
        file containing extra information of river arcs, which is used to identify the river arcs to be dredged.
        You should have this file under the RiverMapper output directory;
        if not, configure RiverMapper to output this file by setting i_DiagnosticOutput
    - region_gdf_file: str, path to the shapefile defining the region of interest,
        in which the dredging is performed. It must have a coordinate reference system (CRS) defined.
    - exclude_region_gdf_file_list: list of str, list of paths to the shapefiles
        defining the regions to be excluded from dredging, e.g., New England (GoME),
        in which auto arcs are used to represent underwater channels, not watershed rivers.
    - output_dir: str, directory to save the dredged mesh and diagnostic files
    '''
    # check inputs:
    if river_extra_info_map_file is None or not os.path.isfile(river_extra_info_map_file):
        raise ValueError('river_extra_info_map_file is required and must be a valid file path.')
    if region_gdf_file is None or not os.path.isfile(region_gdf_file):
        raise ValueError('region_gdf_file is required and must be a valid file path.')
    if output_dir is None:
        raise ValueError('output_dir is required and must be a valid directory path.')
    if exclude_region_gdf_file_list is not None:
        for exclude_region_gdf_file in exclude_region_gdf_file_list:
            if not os.path.isfile(exclude_region_gdf_file):
                raise ValueError(f'exclude_region_gdf_file {exclude_region_gdf_file} must be a valid file path.')

    # Load extra information from the river arcs
    rivers = Rivers(SMS_MAP(river_extra_info_map_file))  # default crs is 'epsg:4326', which is also the default for RiverMapper

    # Define region of interest
    watershed_origional = gpd.read_file(region_gdf_file)
    # Exclude regions
    if exclude_region_gdf_file_list is not None and exclude_region_gdf_file_list != []:
        watershed = deepcopy(watershed_origional)
        for exclude_region_gdf_file in exclude_region_gdf_file_list:
            watershed = gpd.overlay(
                watershed,
                gpd.read_file(exclude_region_gdf_file).to_crs(watershed_origional.crs),
                how='difference'
            )

    # Dredge the river transects
    hgrid_dredged = dredge_river_transects(
        rivers, region_gdf=watershed, hgrid_obj=hgrid_obj,
        min_channel_depth=1.0, output_dir=output_dir,
        measured_from_high_bank=measured_from_high_bank
    )

    return hgrid_dredged
