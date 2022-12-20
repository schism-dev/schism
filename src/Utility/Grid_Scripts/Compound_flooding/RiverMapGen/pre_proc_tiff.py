# %%
import os
from osgeo import gdal
from schism_py_pre_post.Rivermap.make_river_map import Tif2XYZ
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import shapefile
import numpy as np
from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
from pylib import schism_grid
import glob
import pickle


def make_tile_based_on_template(
    newtiles_topleft:np.ndarray,
    template_tile:str,
    background_tile_list:list,
    outdir:str,
):
    # get background tiles' boxes
    bg_box_list = []
    in_bg_tile = np.zeros((len(newtiles_topleft), len(background_tile_list)), dtype=bool)
    for k, tile in enumerate(background_tile_list):
        tif = Tif2XYZ(tile)
        bg_box_list.append([tif.x.min, tif.y.min, tif.x.max, tif.y.max])
        in_bg_tile[:, k] = (newtiles_topleft[:, 0] >= tif.x.min()) * \
                           (newtiles_topleft[:, 0]+0.25 <= tif.x.max()) * \
                           (newtiles_topleft[:, 1]-0.25 >= tif.y.min()) * \
                           (newtiles_topleft[:, 1] <= tif.y.max())

    # read template tile
    tplt_tif = Tif2XYZ(template_tile)
    ds = gdal.Open(template_tile, gdal.GA_Update)
    gt = ds.GetGeoTransform()

    # read background tiles and build interpolators
    f_cache_fname = '/sciclone/schism10/feiye/Cache/f_interp2d.pkl'
    if os.path.exists(f_cache_fname):
        with open(f_cache_fname, 'rb') as file:
            f = pickle.load(file)
    else:
        f = []
        for bg_tif_file in background_tile_list:
            bg_tif = Tif2XYZ(bg_tif_file)
            elev = bg_tif.elev
            elev[np.isnan(elev)] = -9999
            f.append(interpolate.interp2d(bg_tif.x, bg_tif.y, bg_tif.elev, kind='linear'))
        with open(f_cache_fname, 'wb') as file:
            pickle.dump(f, file)

    # loop all additional tiles
    for k, topleft in enumerate(newtiles_topleft):
        print(f'processing tile {k+1} of {len(newtiles_topleft)}')

        lat_name = "{:.2f}".format(abs(topleft[1])).replace('.', 'x')
        lon_name = "{:.2f}".format(abs(topleft[0])).replace('.', 'x')
        tif_out = f'{outdir}/ncei19_n{lat_name}_w0{lon_name}_FromCRM.tif'

        print(f'output file: {tif_out}')
        
        # write new origin
        shift_x = topleft[0] - gt[0]
        shift_y = topleft[1] - gt[3]
        new_gt = list(gt)
        new_gt[0] = topleft[0]
        new_gt[1] = 0.25/(ds.RasterXSize-1)
        new_gt[3] = topleft[1]
        new_gt[5] = -0.25/(ds.RasterYSize-1)
        new_gt = tuple(new_gt)

        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()

        # overwrite elevation with interpolated values from CRM
        i_arr = np.zeros(arr.shape, dtype=bool)
        new_arr = np.ones(arr.shape, dtype='float32') * -9999
        for bg_idx in np.argwhere(in_bg_tile[k]):
            idx = int(bg_idx)
            bg_tif = Tif2XYZ(background_tile_list[idx])
            elev = bg_tif.elev
            elev[np.isnan(elev)] = -9999
            this_f = f[idx]
            tmp = this_f(tplt_tif.x+shift_x, tplt_tif.y+shift_y).astype('float32')
            new_arr[~i_arr] = tmp[~i_arr]
            i_arr += (new_arr > -9999)

        # write new tile
        [rows, cols] = arr.shape
        driver = gdal.GetDriverByName("GTiff")
        outdata = driver.Create(tif_out, cols, rows, 1, gdal.GDT_Float32)
        outdata.SetGeoTransform(new_gt)##sets same geotransform as input
        outdata.SetProjection(ds.GetProjection())##sets same projection as input
        outdata.GetRasterBand(1).WriteArray(np.flipud(new_arr))
        # outdata.GetRasterBand(1).SetNoDataValue(10000.0)##if you want these values transparent
        outdata.FlushCache() ##saves to disk!!
        outdata = None
        band=None

    ds=None

def find_missing_CuDEM_tiles(
    tileindex_file='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/tileindex_NCEI_ninth_Topobathy_2014.shp',
    points=None
):
    # read tile index from CuDEM
    sf = shapefile.Reader(tileindex_file)
    shapes = sf.shapes()

    # CuDEM tile is 0.25 degrees x 0.25 degrees
    pts2parent_tiles = np.c_[np.floor(points[:, 0]*4.0)/4.0, np.ceil(points[:, 1]*4.0)/4.0].view(np.complex128)

    all_needed_tiles = [[x.real, x.imag] for x in np.unique(pts2parent_tiles)]
    available_tiles = [[round(x.bbox[0],2), round(x.bbox[1],2)] for x in shapes]
    needed_tiles = [x for x in all_needed_tiles if x not in available_tiles]

    # np.savetxt('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/all_needed_tiles.xyz', np.array(all_needed_tiles))
    # np.savetxt('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/available_tiles.xyz', np.array(available_tiles))
    # np.savetxt('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/needed_tiles.xyz', np.array(needed_tiles))

    return needed_tiles

    
if __name__ == "__main__":
    tileindex_file = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/tileindex_NCEI_ninth_Topobathy_2014.shp'
    gd_fname = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/I24/hgrid.ll'
    coastal_region_shp = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/coastal_region.shp'
    template_tile_file = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CuDEM/Lonlat/ncei19_n32x25_w081x25_2019v1.tif'
    crm_folder = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CRM/Lonlat/'
    outdir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CuDEM/Lonlat/'

    # read hgrid and find coastal points
    gd = schism_grid(gd_fname)
    coastal = find_node_in_shpfiles(shapefile_names=[coastal_region_shp], gd=gd)
    points = np.c_[gd.x[coastal], gd.y[coastal]]

    # find top-left coordinates of missing CuDEM tiles
    needed_tiles_topleft = find_missing_CuDEM_tiles(tileindex_file=tileindex_file, points=points)

    # make new tiles
    make_tile_based_on_template(
        newtiles_topleft=np.array(needed_tiles_topleft),
        template_tile=template_tile_file,
        background_tile_list=glob.glob(f'{crm_folder}/*.tif'),
        outdir = outdir
    )
