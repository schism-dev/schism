from pylib_essentials.schism_file import read_schism_hgrid_cached, grd2sms
from pathlib import Path
import geopandas as gpd
import numpy as np

wdir = Path('/sciclone/schism10/feiye/Test/RUN02b_JZ/Dredge/')
dredge_depth = 5

hg = read_schism_hgrid_cached(f'{wdir}/hgrid_chart_NCF_loaded.gr3')

# load river polygons
river_polys = gpd.read_file(wdir / 'total_river_polys_clipped_dissolved.shp')
# shrink the polygons by 1 m to exclude bank nodes
river_polys.geometry = river_polys.geometry.to_crs('esri:102008').buffer(-1).to_crs('epsg:4326')

# determine in-channel nodes
hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(hg.x, hg.y), crs='epsg:4326')
joined_gdf = gpd.sjoin(hg_points, river_polys, how="inner", op='within')
idx = joined_gdf.index.to_numpy() 

in_channel = np.zeros_like(hg.dp, dtype=bool)
in_channel[idx] = True
hg.plot(value=in_channel.astype(int), fmt=1)

# dredge the in-channel nodes
hg.dp[idx] += dredge_depth

grd2sms(hg, f'{wdir}/hgrid_chart_NCF_loaded_dredged_{dredge_depth}m.2dm')
hg.save(f'{wdir}/hgrid_chart_NCF_loaded_dredged_{dredge_depth}m.gr3', fmt=1)

pass