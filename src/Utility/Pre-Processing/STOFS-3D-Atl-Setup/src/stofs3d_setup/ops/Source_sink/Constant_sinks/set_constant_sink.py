#!/usr/bin/env python3
import numpy as np
from copy import deepcopy
from pathlib import Path

import geopandas as gpd

from pylib_experimental.schism_file import schism_grid
from pylib_experimental.schism_file import source_sink, TimeHistory


def find_ele_in_shpfile(shapefile_name, hgrid: schism_grid, type=True):
    '''
    Find element index within polygons of a list of shapefiles

    returns:
      inside_ele_mask: a boolean array indicating which elements are inside the polygons
      element_types: a list of strings indicating the type of each element
        if type is true, the element_types assumes the "type" attribute in the parent polygon
        if type is false, the element_types is None
    '''
    polygon_gdf = gpd.read_file(shapefile_name)
    if 'type' not in polygon_gdf.columns:
        # backward compatible (New Orleans used to be the only polygon type)
        # and additional urban polygons used to be assumed to have the same pumping capacity
        polygon_gdf['type'] = 'NewOrleans'

        print(f'Warning: "type" attribute not found in {shapefile_name}, set to "NewOrleans" for all polygons')

    hgrid.compute_ctr()
    pts_gdf = gpd.GeoDataFrame(
        {'ele_idx': np.arange(hgrid.ne)},
        geometry=gpd.points_from_xy(hgrid.xctr, hgrid.yctr),
        crs='esri:102008'  # Albers Equal Area
    )

    joined = gpd.sjoin(pts_gdf, polygon_gdf[["geometry", "type"]], how='left', predicate='within')

    types_by_point = joined.groupby(joined.index)["type"].first()
    # Re-align to ALL points in original order (fill NaN for points outside any polygon)
    ele_types = types_by_point.reindex(pts_gdf.index).to_numpy()

    return ele_types


def set_constant_sink(wdir='./', shapefile_name='levee_4_pump_polys.shp', hgrid: schism_grid = None):
    """
    Set constant sink on land and within polygons defined in shapefile_name
    The sink rate within polygons is set to 0.5 inch/hour, as an estimation of pump capacity
    The sink rate on land is set to 0.5 cm/hour, as an estimation of background infiltration
    """
    def cm_per_hour(cm_h):
        """convert cm/h to m/s"""
        return (cm_h / 100.0) / 3600.0

    def inch_per_hour(inch_h):
        """convert inch/h to m/s"""
        return (inch_h * 0.0254) / 3600.0

    constant_sinks_dict = {
        "water": cm_per_hour(0.02),  # 0.02 cm/h to balance climatological precipitation (~65 in/year in Louisiana)
        "NewOrleans": inch_per_hour(0.5),  # 0.5 in/h ≈ 1.27 cm/h pump capacity
        "urban": cm_per_hour(1.0),  # 1.0 cm/h estimated pump capacity
        "rural": cm_per_hour(0.5),  # 0.5 cm/h infiltration
    }

    # copy datafiles
    gd = deepcopy(hgrid)
    # shapefile_crs = shapefile.Reader(f'{wdir}/{shapefile_name}').crs
    gd.proj(prj0='epsg:4326', prj1='esri:102008')  # project to Albers Equal Area for inside polygon test

    # set pump capacities
    ele_type_array = find_ele_in_shpfile(
        # order matters, latter prevails
        shapefile_name=f'{wdir}/{shapefile_name}',
        hgrid=gd,
    )
    gd.compute_area()

    # diagnostic output
    sink_type_ids = np.zeros((gd.ne, ), dtype=int)  # 0: no sink

    leveed_sinks = np.zeros((gd.ne, ), dtype=float)
    leveed_sink_mask = np.zeros((gd.ne, ), dtype=bool)
    for i, [ele_type, sink_rate] in enumerate(constant_sinks_dict.items()):
        ele_idx = np.where(ele_type_array == ele_type)[0]
        leveed_sinks[ele_idx] = -sink_rate * gd.area[ele_idx]
        leveed_sink_mask[ele_idx] = True
        sink_type_ids[ele_idx] = i + 1  # 0: no sink, 1: water, 2: NewOrleans, etc.

    # set background sink on land
    land = gd.dpe < -2
    background_sink = np.zeros((gd.ne, ), dtype=float)
    background_sink[land] = -0.5/100/3600 * gd.area[land]  # 0.5 cm/h
    const_sinks = np.minimum(background_sink, leveed_sinks)  # take the larger sink (more negative)

    total_sink_eles = (land + leveed_sink_mask).astype(bool)
    # better be a list, for instantiating added_ss
    # index starts from 1
    total_sink_ele_ids = (np.argwhere(total_sink_eles).reshape(-1, ) + 1).tolist()
    # index starts from 0
    total_sink_ele_idx = np.argwhere(total_sink_eles).reshape(-1, )

    # build vsink.th
    vsink = TimeHistory(
        data_array=np.r_[np.c_[np.array([0.0]), const_sinks[total_sink_eles].reshape(1, -1)],
                         np.c_[np.array([100*365*86400]), const_sinks[total_sink_eles].reshape(1, -1)]],
        columns=total_sink_ele_ids,
    )

    # save sink ele ids and sink values for operational use
    np.savetxt(f'{wdir}/pump_sinks.txt', np.c_[total_sink_ele_ids, const_sinks[total_sink_eles]], fmt='%i %10.5f')
    # save as xyz format
    gd.compute_ctr()
    my_ele_xyz = np.c_[gd.xctr[total_sink_ele_idx], gd.yctr[total_sink_ele_idx], sink_type_ids[total_sink_eles]]
    np.savetxt(f'{wdir}/sinks.xyz', my_ele_xyz)

    # build a source_sink object
    const_source_sink = source_sink(vsource=None, vsink=vsink, msource=None)
    # write source sink files, actually only vsink.th and source_sink.in
    const_source_sink.writer(output_dir=wdir)

    return const_source_sink


if __name__ == "__main__":
    my_ss = source_sink.from_files('/sciclone/schism10/feiye/STOFS3D-v8/I202503/Source_sink/')
    eles = np.array(my_ss.source_eles.astype(str))
    # read a json file
    import json
    json_file = Path('/sciclone/schism10/feiye/STOFS3D-v8/I202503/Source_sink/sources.json')
    with open(json_file, 'r') as f:
        ele2fid = json.load(f)
    eles1 = np.array(list(ele2fid.keys()))

    np.all(np.equal(eles, eles1))

    print('Done')
