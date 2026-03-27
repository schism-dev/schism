"""
Test classes and functions in the source_sink module.
"""
from pylib_experimental.schism_file import source_sink
ss = source_sink.from_files(
    "/sciclone/schism10/feiye/STOFS3D-v7.2/Shadow_Forecast"
)


# read json file
import json
with open('/sciclone/schism10/feiye/TEMP/clip_by_polygon/sources.json') as f:
    sources_dict = json.load(f)

tmp = [fid for fids in sources_dict.values() for fid in fids]
if len(tmp) != len(set(tmp)):
    raise ValueError('Duplicated fids in new2fid')

    # print the duplicates
    for fid in set(tmp):
        if tmp.count(fid) > 1:
            print(f'Duplicated fid: {fid}')


import geopandas as gpd
from shapely import get_coordinates
from pylib_experimental.schism_file import source_sink
# from pylib_essentials.schism_file import TimeHistory, source_sink
from schism_py_pre_post.Utilities.import_util import get_hgrid_reader

read_hgrid = get_hgrid_reader()

rundir = '/sciclone/schism10/feiye/STOFS3D-v8/R15c_v7'

hgrid = read_hgrid(f'{rundir}/hgrid.gr3')

# read original source/sink
original_ss = source_sink.from_ncfile(f'{rundir}/source.nc')

region = gpd.read_file('/sciclone/schism10/feiye/TEMP/clip_by_polygon/hires/hires.shp').to_crs("epsg:4326")
region_coords = [get_coordinates(p) for p in region.explode(index_parts=True).exterior]

# split source/sink into inside and outside region
inside_ss, outside_ss = original_ss.clip_by_polygons(hgrid=hgrid, polygons_xy=region_coords,)
print(inside_ss)
print(outside_ss)

print('Done')
