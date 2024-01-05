from pylib_essentials.schism_file import read_schism_hgrid_cached, grd2sms, schism_grid
from pathlib import Path
import geopandas as gpd
import numpy as np

# ----------- inputs -------------
wdir = Path('/sciclone/schism10/feiye/SECOFS/Inputs/I03a/Bathy_edit/Dredge/')
hg_file = Path(f'{wdir}/hgrid.gr3')
dredge_polygon_file = Path(f'{wdir}/total_river_polys_clipped_test.shp')
dredge_depth = 2
# ---------------------------------

hg = schism_grid(str(hg_file))  # epsg:4326

# load river polygons
river_polys = gpd.read_file(dredge_polygon_file)  # epsg:4326, this is from clip_autoarcs.py
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

grd2sms(hg, f'{wdir}/{hg_file.stem}_dredged_{dredge_depth}m.2dm')
hg.save(f'{wdir}/{hg_file.stem}_dredged_{dredge_depth}m.gr3', fmt=1)

pass