import os
import numpy as np
from copy import deepcopy
import geopandas as gpd
import pandas as pd
from RiverMapper.rivers import Rivers
from RiverMapper.SMS import SMS_MAP
from sklearn.neighbors import KDTree
from stofs3d_setup.utils.projection import project_geodataframe


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

    print(f'dredging river transects to maintain a minimum channel depth of {min_channel_depth} m ...')
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
    region_gdf_file_list=None, exclude_region_gdf_file_list=None,
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
    - region_gdf_file_list: list of str, paths to files defining the regions of
        interest in which dredging is performed. Each file must have a coordinate
        reference system (CRS) defined. The regions are combined before dredging.
    - exclude_region_gdf_file_list: list of str, list of paths to the shapefiles
        defining the regions to be excluded from dredging, e.g., New England (GoME),
        in which auto arcs are used to represent underwater channels, not watershed rivers.
    - output_dir: str, directory to save the dredged mesh and diagnostic files
    '''
    # check inputs:
    if river_extra_info_map_file is None or not os.path.isfile(river_extra_info_map_file):
        raise ValueError('river_extra_info_map_file is required and must be a valid file path.')
    if not region_gdf_file_list:
        raise ValueError('region_gdf_file_list is required and cannot be empty.')
    if isinstance(region_gdf_file_list, (str, os.PathLike)):
        raise ValueError('region_gdf_file_list must be a list of file paths.')
    for region_gdf_file in region_gdf_file_list:
        if not os.path.isfile(region_gdf_file):
            raise ValueError(
                f'region_gdf_file {region_gdf_file} must be a valid file path.'
            )
    if output_dir is None:
        raise ValueError('output_dir is required and must be a valid directory path.')
    if exclude_region_gdf_file_list is not None:
        for exclude_region_gdf_file in exclude_region_gdf_file_list:
            if not os.path.isfile(exclude_region_gdf_file):
                raise ValueError(f'exclude_region_gdf_file {exclude_region_gdf_file} must be a valid file path.')

    # Load extra information from the river arcs
    rivers = Rivers(SMS_MAP(river_extra_info_map_file))  # default crs is 'epsg:4326', also the default for RiverMapper

    # Define the combined region of interest. Files are loaded in the supplied
    # order and projected to the CRS of the first file before being combined.
    region_gdf_list = []
    target_crs = None
    for region_gdf_file in region_gdf_file_list:
        region_gdf = gpd.read_file(region_gdf_file)
        if region_gdf.empty:
            raise ValueError(f'region_gdf_file {region_gdf_file} contains no features.')
        if region_gdf.crs is None:
            raise ValueError(f'region_gdf_file {region_gdf_file} has no CRS.')

        if target_crs is None:
            target_crs = region_gdf.crs
        else:
            region_gdf = project_geodataframe(
                region_gdf,
                target_crs,
                f'channel-connectivity inclusion-region projection for {region_gdf_file}',
            )
        region_gdf_list.append(region_gdf[['geometry']])

    watershed = gpd.GeoDataFrame(
        pd.concat(region_gdf_list, ignore_index=True),
        geometry='geometry',
        crs=target_crs,
    ).dissolve().reset_index(drop=True)

    # Exclude regions
    if exclude_region_gdf_file_list is not None and exclude_region_gdf_file_list != []:
        for exclude_region_gdf_file in exclude_region_gdf_file_list:
            watershed = gpd.overlay(
                watershed,
                project_geodataframe(
                    gpd.read_file(exclude_region_gdf_file),
                    watershed.crs,
                    "channel-connectivity exclusion-region projection",
                ),
                how='difference'
            )

    # Dredge the river transects
    hgrid_dredged = dredge_river_transects(
        rivers, region_gdf=watershed, hgrid_obj=hgrid_obj,
        min_channel_depth=min_channel_depth, output_dir=output_dir,
        measured_from_high_bank=measured_from_high_bank
    )

    return hgrid_dredged
