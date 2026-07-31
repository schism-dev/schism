import numpy as np
import os
import shutil
from pathlib import Path
import geopandas as gpd

# pip install git+https://github.com/feiye-vims/schism_py_pre_post.git
from schism_py_pre_post.Shared_modules.hotstart_proc import Hotstart
from schism_py_pre_post.Shared_modules.gen_subregion_ic2 import gen_subregion_ic_stofs3d

# pip install git+https://github.com/wzhengui/pylibs.git
from pylib import schism_grid
from .download_cbp_hotstart import get_cbp_obs_for_stofs3d
from .download_usgs_hotstart import get_usgs_obs_for_stofs3d
from .replace_eta2_aviso import interp_to_points_2d, transform_ll_to_cpp


FALLBACK_OBS_DATE = '2017-12-01'
REQUIRED_OBS_PREFIXES = ('mean_sal_xyz', 'mean_tem_xyz')


def _count_data_lines(path):
    if not path.exists():
        return 0
    with path.open('r', encoding='utf-8') as fin:
        return sum(1 for line in fin if line.strip())


def _obs_file(obs_dir, prefix, date_str):
    return Path(obs_dir) / f'{prefix}_{date_str}'


def _obs_files_are_ready(obs_dir, date_str, fallback_obs_dir, min_fallback_fraction=0.25):
    for prefix in REQUIRED_OBS_PREFIXES:
        obs_file = _obs_file(obs_dir, prefix, date_str)
        fallback_file = _obs_file(fallback_obs_dir, prefix, FALLBACK_OBS_DATE)
        fallback_lines = _count_data_lines(fallback_file)
        min_lines = max(1, int(fallback_lines * min_fallback_fraction)) if fallback_lines > 0 else 1
        if _count_data_lines(obs_file) < min_lines:
            return False
    return True


def _use_fallback_obs(obs_dir, date_str, fallback_obs_dir):
    obs_dir = Path(obs_dir)
    obs_dir.mkdir(parents=True, exist_ok=True)

    print('\n' + '!' * 78)
    print('WARNING: LIVE USGS OBS DOWNLOAD FAILED OR PRODUCED INCOMPLETE FILES.')
    print(f'WARNING: USING FALLBACK HOTSTART OBS FROM {fallback_obs_dir}')
    print(f'WARNING: RENAMING FALLBACK OBS DATE {FALLBACK_OBS_DATE} TO RUN DATE {date_str}.')
    print('!' * 78 + '\n')

    for prefix in REQUIRED_OBS_PREFIXES:
        source = _obs_file(fallback_obs_dir, prefix, FALLBACK_OBS_DATE)
        destination = _obs_file(obs_dir, prefix, date_str)
        if not source.exists():
            raise FileNotFoundError(f'Fallback obs file not found: {source}')
        if destination.exists() or destination.is_symlink():
            destination.unlink()
        shutil.copy2(source, destination)


def _read_obs_points(obs_file):
    obs_file = Path(obs_file)
    rows = []
    if not obs_file.exists():
        return rows

    with obs_file.open('r', encoding='utf-8') as fin:
        for line in fin:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                lon = float(parts[0])
                lat = float(parts[1])
                value = float(parts[2])
            except ValueError:
                continue
            if not np.isfinite(lon) or not np.isfinite(lat) or not np.isfinite(value):
                continue
            rows.append((lon, lat, value, parts[3]))
    return rows


def _obs_source(station_id):
    return 'USGS' if str(station_id).isdigit() else 'CBP'


def plot_hotstart_obs(obs_dir, date_str, outfilename=None):
    obs_dir = Path(obs_dir)
    outfilename = Path(outfilename or obs_dir / f'obs_scatter_salinity_temperature_{date_str}.png')
    datasets = [
        ('Temperature', 'deg C', _read_obs_points(_obs_file(obs_dir, 'mean_tem_xyz', date_str))),
        ('Salinity', 'psu', _read_obs_points(_obs_file(obs_dir, 'mean_sal_xyz', date_str))),
    ]
    all_rows = [row for _, _, rows in datasets for row in rows]
    if not all_rows:
        print('WARNING: no observations found for diagnostic scatter plot')
        return None

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    lon = np.array([row[0] for row in all_rows])
    lat = np.array([row[1] for row in all_rows])
    lon_pad = max(1.0, 0.08 * (lon.max() - lon.min()))
    lat_pad = max(1.0, 0.08 * (lat.max() - lat.min()))
    extent = [lon.min() - lon_pad, lon.max() + lon_pad, lat.min() - lat_pad, lat.max() + lat_pad]

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        projection = ccrs.PlateCarree()
        fig, axes = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': projection}, constrained_layout=True)
        for ax in axes:
            ax.set_extent(extent, crs=projection)
            ax.add_feature(cfeature.LAND, facecolor='0.86', edgecolor='none', zorder=0)
            ax.add_feature(cfeature.OCEAN, facecolor='white', edgecolor='none', zorder=0)
            ax.add_feature(cfeature.COASTLINE, edgecolor='0.35', linewidth=0.7, zorder=1)
            ax.add_feature(cfeature.BORDERS, edgecolor='0.60', linewidth=0.4, zorder=1)
            ax.gridlines(draw_labels=True, linewidth=0.3, color='0.70', alpha=0.7)
        transform = projection
    except Exception as exc:
        print(f'WARNING: cartopy basemap failed ({exc}); writing lon/lat scatter without shoreline')
        fig, axes = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True)
        for ax in axes:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
            ax.set_facecolor('0.94')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            ax.grid(True, linewidth=0.3, color='0.80')
        transform = None

    markers = {'USGS': 'o', 'CBP': '^'}
    for ax, (title, units, rows) in zip(axes, datasets):
        ax.set_title(f'{title} observations, {date_str}')
        if not rows:
            ax.text(0.5, 0.5, 'No observations', ha='center', va='center', transform=ax.transAxes)
            continue
        values_all = [row[2] for row in rows]
        vmin = min(values_all)
        vmax = max(values_all)
        scatter_for_colorbar = None
        for source, marker in markers.items():
            source_rows = [row for row in rows if _obs_source(row[3]) == source]
            if not source_rows:
                continue
            x = [row[0] for row in source_rows]
            y = [row[1] for row in source_rows]
            values = [row[2] for row in source_rows]
            kwargs = {
                'c': values,
                's': 26,
                'marker': marker,
                'edgecolors': 'black',
                'linewidths': 0.25,
                'label': f'{source} ({len(source_rows)})',
                'vmin': vmin,
                'vmax': vmax,
                'zorder': 3,
            }
            if transform is not None:
                kwargs['transform'] = transform
            scatter_for_colorbar = ax.scatter(x, y, **kwargs)
        if scatter_for_colorbar is not None:
            cbar = fig.colorbar(scatter_for_colorbar, ax=ax, shrink=0.78)
            cbar.set_label(units)
        ax.legend(loc='lower left', fontsize=9, frameon=True)

    outfilename.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfilename, dpi=180)
    plt.close(fig)
    print(f'wrote hotstart observation diagnostic plot: {outfilename}')
    return outfilename


def find_nodes_in_shapefile(hgrid, shapefile_fname, crs):
    """
    Find nodes in polygons defined in a shapefile
    """
    hg_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(hgrid.x, hgrid.y), crs=crs)
    polygon_gdf = gpd.read_file(shapefile_fname)
    if 'type' in polygon_gdf.columns:
        polygon_gdf = polygon_gdf[polygon_gdf['type'].fillna('').str.lower() != 'water'].copy()
    if hg_points.crs != polygon_gdf.crs:
        hg_points = hg_points.to_crs(polygon_gdf.crs)
    ids = gpd.sjoin(hg_points, polygon_gdf, predicate='within', how='inner').index
    idx_mask = np.zeros(hgrid.dp.shape, dtype=bool)
    idx_mask[ids] = True

    return ids, idx_mask


def gen_elev_ic(hgrid=None, h0=0.1, city_shape_fnames=None, aviso_file=None):
    '''
    set initial elevation: 0 in the ocean, h0 below ground on higher grounds and in cities
    city_shape_fname: the shapefile must be in the same projection as the hgrid
    in this case lon/lat (because aviso data is in lon/lat)
    '''
    in_city = np.zeros_like(hgrid.dp, dtype=bool)
    for shp in city_shape_fnames:
        _, idx_mask = find_nodes_in_shapefile(hgrid, shp, crs='EPSG:4326')
        in_city = np.logical_or(in_city, idx_mask)

    # in_city = find_points_in_polyshp(pt_xy=np.c_[hgrid.x, hgrid.y], shapefile_names=city_shape_fnames)

    hgrid.save('incity.gr3', value=in_city.copy().astype(int))

    high_ground = hgrid.dp < 0

    elev_ic = np.zeros(hgrid.dp.shape, dtype=float)

    # set coastal points
    land = np.logical_or(high_ground, in_city)
    elev_ic[land] = - hgrid.dp[land] - h0

    # set ocean points
    import xarray as xr
    ocean = np.logical_not(land)
    x1, y1 = transform_ll_to_cpp(hgrid.x, hgrid.y)

    if aviso_file is None:
        print('no aviso file provided, setting elev to 0 in the ocean')
    else:
        print(f'reading aviso data from {aviso_file} for ocean elev')
        # read aviso data
        ds = xr.open_dataset(aviso_file)
        lon2 = ds.longitude.values
        lat2 = ds.latitude.values
        x2, y2 = transform_ll_to_cpp(lon2, lat2)
        ssh = np.squeeze(ds.adt.values[0, :, :])
        ds.close()
        # interpolate aviso data onto model grid
        ssh_int = interp_to_points_2d(y2, x2, np.c_[y1, x1], ssh)

        elev_ic[ocean] = ssh_int[ocean] - 0.42  # nEw, uniform shift for STOFS-3D-Atl v7

    return elev_ic


def tweak_stofs3d_hotstart(
    wdir='./', hotstart_date_str='2021-05-01', city_shapefile_names=[], aviso_file='aviso.nc',
    hycom_TS_file='TS_1.nc', hycom_hot_file='hotstart.nc.hycom',
):
    '''
    The "wdir" should contain the following files:
        hgrid.gr3 (should be the same as hgrid.ll for STOFS3D), vgrid.in
        TS_1.nc: Temperature and Salinity from HYCOM, the first time in the nc file should match model start time.
                 This file should already be prepared in the previous step of hotstart.nc generation
        hotstart.nc.hycom: hotstart.nc based on HYCOM only, renamed with the '.hycom' suffix.
        aviso.nc: sea surface height above geoid, in lon/lat

    You may select pre-defined city polygons to force the water level to be below ground;
    only specify the file basename, the script will copy these files to your wdir.

    Use scripts in ../AVISO/ to download aviso.nc
    '''
    griddir = wdir
    output_obs_dir = f'{wdir}/Obs/'  # observations will be downloaded to this directory
    my_hot_file = f'{wdir}/hotstart.nc'  # this is the final hotstart.nc file
    hycom_TS_file = f'{wdir}/TS_1.nc'  # this file should already be prepared in the previous step of hotstart.nc generation from HYCOM
    hycom_hot_file = f'{wdir}/hotstart.nc.hycom'  # this file should already be prepared in the previous step of hotstart.nc generation from HYCOM, renamed with the '.hycom' suffix.
    fallback_obs_dir = Path(os.path.dirname(os.path.abspath(__file__))) / 'Obs_2017_12_01'

    # copy datafiles
    mydir = os.path.dirname(os.path.abspath(__file__))
    for shp in city_shapefile_names:
        shp_basename = Path(shp).stem
        if not os.path.exists(f'{wdir}/{shp_basename}.shp'):
            os.system(f'cp {mydir}/{shp_basename}.* {wdir}')

    if _obs_files_are_ready(output_obs_dir, hotstart_date_str, fallback_obs_dir):
        print('observation files already exist, skip downloading')
    else:
        print(f'downloading obs from USGS to {output_obs_dir}')
        try:
            get_usgs_obs_for_stofs3d(
                outdir=output_obs_dir,
                start_date_str=hotstart_date_str,
                hgrid_path=f'{wdir}/hgrid.gr3',
            )

            print('downloading coastal obs from CBP')
            get_cbp_obs_for_stofs3d(outdir=output_obs_dir, sample_time=hotstart_date_str, varname=['sal', 'tem'])
        except Exception as exc:
            print(f'WARNING: live observation download raised an exception: {exc}')

        if not _obs_files_are_ready(output_obs_dir, hotstart_date_str, fallback_obs_dir):
            _use_fallback_obs(output_obs_dir, hotstart_date_str, fallback_obs_dir)

    plot_hotstart_obs(
        output_obs_dir,
        hotstart_date_str,
        Path(wdir) / f'obs_scatter_salinity_temperature_{hotstart_date_str}.png',
    )

    print('interpolating obs onto model grid')
    gen_subregion_ic_stofs3d(wdir=wdir, obsdir=output_obs_dir, hycom_TS_file=hycom_TS_file, date_str=hotstart_date_str)

    print(f'copying {hycom_TS_file} to {my_hot_file}, which will be modified')
    if os.path.exists(my_hot_file):
        os.system(f"rm {my_hot_file}")
    os.system(f"cp {hycom_hot_file} {my_hot_file}")

    my_hot = Hotstart(
        grid_info=griddir,
        hot_file=my_hot_file
    )

    # increase T by 1 oC
    # my_hot.tr_nd.val[:, :, 0] += 1.0
    # my_hot.tr_nd0.val[:, :, 0] += 1.0
    # my_hot.tr_el.val[:, :, 0] += 1.0

    print('tweaking coastal values based on obs')
    for i, var in enumerate(['tem', 'sal']):
        hg = schism_grid(f'{wdir}/ecgc_coastal_{var}.gr3')  # this file is from get*_obs_for_stofs3d
        idx = hg.dp > -9998
        for k in range(my_hot.grid.vgrid.nvrt):
            my_hot.tr_nd.val[idx, k, i] = hg.dp[idx]

    print('setting salinity to 0 on higher grounds')
    rat = np.maximum(np.minimum(1.0, (my_hot.grid.hgrid.dp + 3.0) / 3.0), 0.0)  # linearly varying from 0 to 3 m
    my_hot.tr_nd.val[:, :, 1] *= np.transpose(np.tile(rat, (my_hot.grid.vgrid.nvrt, 1)))
    my_hot.trnd_propogate()  # propogate trnd values to trnd0 and tr_el

    print('setting initial elevation: aviso values in the ocean, just below ground on higher grounds and in cities')
    my_hot.eta2.val[:] = gen_elev_ic(
        hgrid=my_hot.grid.hgrid, h0=0.1,
        city_shape_fnames=[f'{wdir}/{x}' for x in city_shapefile_names],
        aviso_file=aviso_file
    )

    print('writing the final hotstart.nc')
    my_hot.writer(fname=my_hot_file)
    return True

if __name__ == "__main__":
    # Modify hycom-based hotstart.nc.hycom with coastal observation values of T and S, and a better initial elevation
    tweak_stofs3d_hotstart(
        wdir='./',
        hotstart_date_str='2024-03-05',
        city_shapefile_names=["levee_pump_polys_2026_with_poly_type_lonlat.shp"],
        aviso_file='aviso.nc'
    )
