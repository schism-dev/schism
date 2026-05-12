#!/usr/bin/env python

"""This script is to adjust the vsource.th file based on USGS data."""

import os

os.environ["MPLBACKEND"] = "Agg"   # or set in your shell: export MPLBACKEND=Agg
os.environ.pop("DISPLAY", None)    # optional: avoid accidental GUI use

from pathlib import Path
import time
from datetime import datetime  # , timedelta
from glob import glob
import json
import pickle

from scipy.spatial import cKDTree
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import xarray as xr
import geopandas as gpd

from schism_py_pre_post.Grid.Prop import Prop
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Utilities.util import b_in_a
from schism_py_pre_post.Download.download_usgs import \
    download_stations, usgs_var_dict, convert_to_ObsData, \
    get_usgs_stations_from_state, detect_data_gap
from stofs3d_setup.utils.utils import STOFS3D_ATL_STATES
from pylib import schism_grid


cache_folder = '/sciclone/schism10/feiye/Cache/'


def read_nwm_data(
    data_dir=None, fids=None,
    date_range=pd.date_range("2015-06-01 00:00:00", "2015-08-01 00:00:00", freq="3H")
):
    '''read NWM data for a list of fids and a date range'''
    if fids is None:
        fids = []

    # assemble data file list
    datafiles = []
    for this_date in date_range:
        files = glob(f'{data_dir}/*{datetime.strftime(this_date, "%Y%m%d%H")}*')
        if len(files) == 0:
            raise FileNotFoundError(f'no file found for {this_date}')
        if len(files) > 1:
            raise ValueError(f'multiple files found for {this_date}')
        datafiles.append(files[0])

    # read data
    stream_flow = np.empty((len(datafiles), len(fids)))
    stream_flow[:] = np.nan
    for i, this_file in enumerate(datafiles):
        print(f'reading {this_file}\n')
        with xr.open_dataset(this_file) as ds:
            fid_all = ds.feature_id.data
            fids_in_all = b_in_a(a=fid_all, b=fids)
            if np.ma.is_masked(fids_in_all):
                raise ValueError('failed to find some segments in NWM data')
            stream_flow[i, :] = ds.streamflow.data[fids_in_all]

    return stream_flow


class UsgsSt():
    '''simple class to store USGS station ID and the nearby NWM segment ID.'''
    def __init__(self, st_id=None, nearby_nwm_fid=None):
        self.st_id = st_id
        self.nearby_nwm_fid = nearby_nwm_fid


class Vsource:
    '''simple class to store vsource information.'''
    def __init__(self, xyz=None, hgrid_ie=0, nwm_fid=0, df=pd.DataFrame()):
        if xyz is None:
            xyz = [0.0, 0.0, 0.0]

        self.xyz = xyz
        self.usgs_st = []
        self.nwm_fid = nwm_fid

        self.upstream_path = SearchPath()
        self.downstream_path = SearchPath()
        self.hgrid_ie = hgrid_ie
        self.df = df  # time series from vsource.th

    @property
    def df(self):
        '''getter for the time series data frame of the vsource.'''
        return self._df

    @df.setter
    def df(self, dataframe):
        if isinstance(dataframe, pd.DataFrame):
            if 'datetime' in dataframe.keys():
                dataframe = dataframe.set_index(
                    pd.DatetimeIndex(pd.to_datetime(dataframe["datetime"], utc=True)))
                self._df = dataframe
        else:
            raise ValueError("needs an instance of the 'TimeHistory' class.")


class SearchPath:
    '''
    This class facilitates the search along the upstream or downstream path of a vsource.
    '''
    def __init__(self):
        self._seg_id = []
        self._order = []

    def forward_path(self, seg_id, seg_order):
        '''forward path by one segment.'''
        self._seg_id.append(seg_id)
        self._order.append(seg_order)

    def reset_path(self):
        '''reset the path to empty'''
        self._seg_id = []
        self._order = []

    def writer(self, fname, shp_gpd=None):
        '''write the path to a file.'''
        with open(fname, 'w', encoding='utf-8') as fout:
            fout.write('lon lat \n')
            for fid in self._seg_id:
                fout.write(str(shp_gpd.loc[fid].lon) + " " + str(shp_gpd.loc[fid].lat) + "\n")


def find_usgs_along_nwm(
    iupstream=True, starting_seg_id=None, vsource=None, total_seg_length=0.0, order=0, nwm_shp=None
):
    """
    This is a recursive function to find the up/down-stream USGS stations,
    The search stops when the total length of the path exceeds 30 km.
    The order of the stream is also recorded,
    but not used as a stopping criterion at this moment.
    """
    SEARCH_DIST = 30e3  # search distance limit

    if starting_seg_id is None:
        return  # some starting segs may not be in nwm_shp, which is subsetted

    if starting_seg_id not in nwm_shp.index:
        print(f'starting segment {starting_seg_id} not found in the shapefile, skipping...')
        return

    if iupstream:
        vsource.upstream_path.forward_path(seg_id=starting_seg_id, seg_order=order)
    else:
        vsource.downstream_path.forward_path(seg_id=starting_seg_id, seg_order=order)

    total_seg_length += nwm_shp.loc[starting_seg_id].Length
    if total_seg_length > SEARCH_DIST:  # stop search when the path exceeds 30 km
        return

    if nwm_shp.loc[starting_seg_id].gages is not None and nwm_shp.loc[starting_seg_id].gages != '':  # gage found
        st = UsgsSt(st_id=nwm_shp.loc[starting_seg_id].gages, nearby_nwm_fid=starting_seg_id)
        vsource.usgs_st.append(st)
        return
    else:  # search upstream or downstream
        if iupstream:
            uids = nwm_shp.loc[starting_seg_id]['from']
        else:
            uids = [nwm_shp.loc[starting_seg_id]['to']]  # 'to' is a single value, make it a list for the for loop

        for uid in uids:
            find_usgs_along_nwm(iupstream=iupstream, starting_seg_id=uid,
                                vsource=vsource, total_seg_length=total_seg_length, order=order+1, nwm_shp=nwm_shp)


def prepare_usgs_stations(states=None, diag_output=None):
    '''
    Prepare USGS stations for a spatial range larger than the model domain,
    so that stations within the search range can be selected for each vsource.

    states: list of states to search for USGS stations, e.g., ['LA', 'MS', 'TX']
    '''
    if states is None:
        raise ValueError('a list of states must be provided')

    # get station info
    station_info = get_usgs_stations_from_state(
        states=states, parameter=usgs_var_dict["streamflow"]["id"])

    station_ids = station_info["site_no"].to_numpy()
    station_lon = station_info['dec_long_va'].to_numpy()
    station_lat = station_info['dec_lat_va'].to_numpy()

    # get rid of invalid stations
    valid = (station_lon != '') & (station_lat != '')
    station_ids = station_ids[valid]
    station_lon = station_lon[valid].astype(float)
    station_lat = station_lat[valid].astype(float)

    # remove duplicates
    _, unique_idx = np.unique(station_ids, return_index=True)
    station_ids = station_ids[unique_idx]
    station_lon = station_lon[unique_idx]
    station_lat = station_lat[unique_idx]

    if diag_output is not None:
        with open(diag_output, 'w', encoding='utf-8') as f:
            f.write('lon,lat,id\n')
            for i, _ in enumerate(station_ids):
                f.write(f'{station_lon[i]},{station_lat[i]},{station_ids[i]}\n')

    return station_ids, np.c_[station_lon, station_lat]


def associate_poi_with_nwm(
    nwm_shp, poi: np.ndarray,
    poi_names: list = None,
    poi_label: str = 'poi',
    diag_output: str = None
) -> gpd.geodataframe.GeoDataFrame:
    '''
    Associate USGS stations with NWM segment fid.
    This is done by finding the nearest NWM segment to each USGS station within a threshold.
    Then, it ensures only one USGS station is associated with each NWM segment.

    Note that some NWM segments have already associated USGS stations in the shapefile
    by the attribute 'gages', which are from NWM v1. In this case, the gages will be
    overwritten with the newly associated USGS stations.

    Input:
    nwm_shp: geopandas dataframe of NWM hydrofabric shapefile
    poi: point of interest, a numpy array of [lon, lat], e.g., of USGS stations
    '''

    if poi_names is None:
        poi_names = [f'{i}' for i in range(len(poi))]
    poi_names = np.array(poi_names)

    # Build KDTree from segment centroids
    seg_centers = np.c_[nwm_shp.geometry.centroid.x, nwm_shp.geometry.centroid.y]
    seg_tree = cKDTree(seg_centers)

    # query k nearest centroid neighbors
    k = 10
    dist, idx = seg_tree.query(poi, k=k)   # (n_poi, k)

    # enforce threshold (10 km ≈ 0.1°)
    dist[dist > 0.1] = np.inf
    idx[dist == np.inf] = -1

    poi_points = gpd.points_from_xy(poi[:, 0], poi[:, 1])

    poi2nwm_idx = np.full(len(poi), -1, dtype=int)
    poi2nwm_distance = np.full(len(poi), np.inf)

    for i in range(len(poi)):
        candidates = []
        for j in range(k):
            seg_pos = idx[i, j]
            if seg_pos == -1:
                continue
            seg_idx = nwm_shp.index[seg_pos]
            # true geometry distance
            geom_dist = nwm_shp.loc[seg_idx, "geometry"].distance(poi_points[i])
            if geom_dist <= 0.01:   # 1 km cutoff
                candidates.append((seg_idx, geom_dist))

        if candidates:
            # pick the closest by geometry distance
            seg_idx, geom_dist = min(candidates, key=lambda x: x[1])
            poi2nwm_idx[i] = nwm_shp.index.get_loc(seg_idx)  # store positional index
            poi2nwm_distance[i] = geom_dist

    # deal with multiple poi's associated with the same nwm segment,
    # retain the one with the smallest distance
    # Build (poi -> seg) matches that passed screening
    valid = np.where(poi2nwm_idx != -1)[0]
    matches = pd.DataFrame({
        "poi": [poi_names[i] for i in valid],
        "seg_pos": [poi2nwm_idx[i] for i in valid],          # position in nwm_shp.index
        "dist": [poi2nwm_distance[i] for i in valid],
    })
    seg_index = nwm_shp.index.to_numpy()
    matches["seg_idx"] = seg_index[matches["seg_pos"].values]

    # --- WARN: multiple POIs on the same segment ---
    seg_dupes = matches[matches.duplicated("seg_idx", keep=False)].sort_values(["seg_idx", "dist"])
    for seg_idx, grp in seg_dupes.groupby("seg_idx"):
        fid = nwm_shp.at[seg_idx, "featureID"] if "featureID" in nwm_shp.columns else seg_idx
        print(f"warning: NWM segment {fid} matched {len(grp)} POIs:")
        for _, r in grp.iterrows():
            print(f"    poi={r.poi}, dist={r.dist:.6f}")
        print(f"    choosing the closest one: poi={grp.iloc[0].poi}, dist={grp.iloc[0].dist:.6f}\n")
    # --- Keep the closest POI per segment ---
    best_seg = matches.sort_values(["seg_idx", "dist"]).groupby("seg_idx", as_index=False).first()
    # reset/assign column
    if poi_label not in nwm_shp.columns:
        nwm_shp[poi_label] = None
    nwm_shp[poi_label] = None
    nwm_shp.loc[best_seg["seg_idx"].values, poi_label] = best_seg["poi"].values

    idx = nwm_shp[poi_label].notnull()
    print(f'{sum(idx)} {poi_label} associated with NWM segments')

    # save diagnostic file
    if diag_output is not None:
        with open(diag_output, 'w', encoding='utf-8') as f:
            for i, this_idx in enumerate(poi2nwm_idx):
                if this_idx != -1:
                    idx = nwm_shp.index[this_idx]
                    f.write(f'{poi[i, 0]}, {poi[i, 1]}, {poi_names[i]}, {nwm_shp.loc[idx].featureID}\n')

    return nwm_shp


def preprocess_nwm_shp(f_shapefile):
    '''
    Preprocess the NWM hydrofabric shapefile.
    Add a "from" property to each segment, i.e., the featureID of the upstream segments.
    '''
    t1 = time.time()

    if os.path.exists(f'{f_shapefile}.pkl'):  # cached
        with open(f'{f_shapefile}.pkl', 'rb') as file:
            nwm_shp = pickle.load(file)
    else:
        nwm_shp = gpd.read_file(f_shapefile)

        # append "from" property to each seg
        nwm_shp['from'] = None

        index_g2l = {}  # global to local index
        for i, rec in nwm_shp.iterrows():
            index_g2l[rec.featureID] = i

        from_list = [[] for _ in range(len(nwm_shp))]  # directly use nwm_shp.loc is slow
        for i, rec in nwm_shp.iterrows():
            try:
                from_list[index_g2l[rec.to]].append(rec.featureID)  # index is featureID
            except KeyError:
                pass  # index_g2l[rec.to] can be None, then 0 is pre-assigned
        nwm_shp['from'] = from_list

        # set featureID as index
        nwm_shp = nwm_shp.set_index('featureID')
        # duplicate featureID as a column
        nwm_shp['featureID'] = nwm_shp.index

        with open(f'{f_shapefile}.pkl', 'wb') as file:
            pickle.dump(nwm_shp, file)

    print(f'---------------loading` and processing shapefile took: {time.time()-t1} s ---------------\n')
    return nwm_shp


def _fingerprint(path: str) -> dict:
    import hashlib
    """Quick file fingerprint to detect source changes."""
    st = os.stat(path)
    return {
        "size": st.st_size,
        "mtime": int(st.st_mtime),
        "sha1_1k": hashlib.sha1(open(path, "rb").read(1024)).hexdigest(),
    }


def preprocess_nwm_shp2(f_shapefile: str) -> gpd.GeoDataFrame:
    """
    Preprocess the NWM hydrofabric shapefile.
    Adds a 'from' property (list of upstream featureIDs) to each segment.
    Caches the processed GeoDataFrame as a parquet file for faster reloads.
    """

    t1 = time.time()
    cache_file = f"{cache_folder}/{Path(f_shapefile).stem}_added_from.parquet"
    cache_meta = f"{cache_folder}/{Path(f_shapefile).stem}_added_from.metadata.json"

    # Compute source fingerprint
    src_fp = _fingerprint(f_shapefile)

    # Load from cache if valid
    if os.path.exists(cache_meta):
        try:
            meta = json.load(open(cache_meta))
            if meta.get("source_fingerprint") == src_fp and os.path.exists(cache_file):
                gdf = gpd.read_parquet(cache_file)
                print(f"Loaded cached parquet ({time.time()-t1:.1f}s)")
                return gdf
        except Exception:
            pass  # ignore broken metadata and rebuild

    # --- rebuild ---
    print("Rebuilding parquet cache ...")
    gdf = gpd.read_file(f_shapefile)

    # Ensure required fields exist
    required = {"featureID", "to"}
    missing = required - set(gdf.columns)
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # Build 'from' relationships vectorized
    upstream_map = gdf.groupby("to")["featureID"].apply(list)
    gdf["from"] = gdf["featureID"].map(upstream_map).apply(
        lambda x: x if isinstance(x, list) else []
    )

    # Set featureID as index (and keep as column)
    gdf = gdf.set_index("featureID", drop=False)

    # Save to GPKG
    gdf.to_parquet(cache_file, index=False)

    # Write metadata
    json.dump(
        {
            "source": os.path.abspath(f_shapefile),
            "source_fingerprint": src_fp,
            "cached_at": int(time.time()),
            "format": "parquet",
            "columns": list(gdf.columns),
            "crs": gdf.crs.to_wkt() if gdf.crs else None,
        },
        open(cache_meta, "w"),
        indent=2,
    )

    print(f"Processed and cached as parquet ({time.time()-t1:.1f}s)")
    return gdf


def test():
    '''test function for debugging'''

    # my_obs = ObsData(usgs_data_cache_fname, from_raw_data=False)

    # ids = [station.id for station in my_obs.stations]

    # with open(f'{original_ss_dir}/ele_fid_dict', encoding='utf-8') as json_file:
    #     mysrc_nwm_fid = json.load(json_file)

    with open('/sciclone/schism10/lcui01/schism20/ICOGS/ICOGS3D/Scripts/RUN24/NWM/sources.json',
              encoding='utf-8') as json_file:
        mysrc_nwm_fid = json.load(json_file)
    n = 0
    for key, value in mysrc_nwm_fid.items():
        if len(value) > 1:
            n += 1
        mysrc_nwm_fid[key] = int(value[0])


CFS_TO_M3S = 0.028316846592


def infer_step(index: pd.DatetimeIndex) -> pd.Timedelta:
    """Infer a representative step for an (almost) regular index."""
    if len(index) < 2:
        return pd.Timedelta(0)
    d = pd.Series(index).diff().dropna()
    # median is robust if a few irregularities exist
    return pd.to_timedelta(d.median())


def blend_usgs_to_vsource(
    vs_time_index: pd.DatetimeIndex,
    vsource_col_m3s: np.ndarray,
    usgs_df: pd.DataFrame,
    *,
    value_col: int = 1,
) -> np.ndarray:
    """
    Replace vsource with USGS instantaneous values where available,
    Reindex to vsource time index, retain any data gap
    """
    if usgs_df is None or usgs_df.empty:
        return vsource_col_m3s.copy()

    vs_idx = pd.DatetimeIndex(vs_time_index)
    usgs_idx = pd.DatetimeIndex(usgs_df.index)

    # --- build USGS series in m^3/s ---
    s_obs = pd.Series(usgs_df.iloc[:, value_col].astype(float).values, index=usgs_idx)
    s_obs = s_obs.replace([np.inf, -np.inf], np.nan) * CFS_TO_M3S
    if s_obs.dropna().empty:
        return vsource_col_m3s.copy()

    # --- exact match only (no nearest, no interpolation) ---
    s_vs = s_obs.reindex(vs_idx)  # NaN where no exact USGS stamp exists

    # --- blend: use USGS where present; otherwise keep original vsource ---
    out = np.where(s_vs.notna().values, s_vs.values, vsource_col_m3s)

    return out


def source_nwm2usgs(
    start_time_str="2017-12-01 00:00:00",
    states=['LA', 'MS', 'TX'],
    f_shapefile="/sciclone/schism10/feiye/schism20/REPO/NWM/Shapefiles/ecgc/ecgc.shp",
    original_ss_dir='/sciclone/schism10/feiye/Requests/RUN02a_JZ/src/NWM/',
    nwm_data_dir='/sciclone/schism10/whuang07/schism20/NWM_v2.1/',
    output_dir='/sciclone/schism10/feiye/Requests/RUN02a_JZ/src/NWM/USGS_adjusted_sources/',
):
    '''
    This script is to adjust the vsource.th file based on USGS data.
    The adjustment is done by comparing the NWM data near the vsource and the USGS data near the vsource.
    The USGS data is downloaded from USGS website using the API.
    The NWM data is pre-downloaded during the original NWM to source/sink conversion by gen_source_sink.py.
    See the Atlas paper draft for details of how USGS stations are selected for each vsource.
    In short, the USGS stations are selected along the upstream and downstream branches of the vsource.
    The NWM shapefile provides the information of the upstream and downstream branches,
    as well as the location of the USGS stations.
    However, only NWM v1's shapefile has the USGS station information.
    This is acceptable because the segment IDs are the same between NWM v1 and later versions.

    The output is adjusted vsource.th written to output_dir.
    '''

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f'{output_dir}/Diag/', exist_ok=True)
    os.makedirs(f'{output_dir}/Svgs/', exist_ok=True)

    cache_name = f'usgs_states_{"_".join(np.sort(states))}.txt'
    if os.path.exists(f'{cache_folder}/{cache_name}'):
        usgs_station_df = pd.read_csv(f'{cache_folder}/{cache_name}', dtype=object)
        usgs_stations = usgs_station_df['id']
        usgs_stations_coords = usgs_station_df[['lon', 'lat']].to_numpy().astype(float)
    else:
        usgs_stations, usgs_stations_coords = prepare_usgs_stations(
            states=states,
            diag_output=f'{output_dir}/usgs_stations.csv',
        )

    nwm_shp = preprocess_nwm_shp(f_shapefile)

    # associate USGS stations with NWM segment fid
    nwm_shp = associate_poi_with_nwm(
        nwm_shp, poi=usgs_stations_coords, poi_names=usgs_stations.tolist(),
        poi_label='gages', diag_output=f'{output_dir}/auto_nearby_gages.txt')

    # manually associate some USGS stations with NWM segments,
    # set the NWM segment to be the one associated with the vsource injection
    manual_nwm2usgs = {
        19406836: '07381490',  # Atchafalaya River at Simmesport, LA
        15708755: '02489500',  # Pearl River at Bogalusa, LA
        18928090: '07375175',  # Bogue Falaya River at Boston St, Covington, LA
    }
    for nwm_featureID, usgs_id in manual_nwm2usgs.items():
        nwm_shp.loc[nwm_shp['featureID'] == nwm_featureID, 'gages'] = usgs_id

    idx = nwm_shp['gages'].notnull()
    print(f'{nwm_shp.loc[idx]}')

    # -----------------------------------------------------------------------------------
    # ------load vsource------
    # -----------------------------------------------------------------------------------

    # read hgrid to determine vsource locations
    t1 = time.time()
    hg = schism_grid(f'{original_ss_dir}/hgrid.gr3')  # hg.save()
    hg.compute_ctr()

    # read featureID at vsource elements
    if os.path.exists(f'{original_ss_dir}/ele_fid_dict.json'):
        with open(f'{original_ss_dir}/ele_fid_dict.json', encoding='utf-8') as json_file:
            mysrc_nwm_fid = json.load(json_file)
    elif os.path.exists(f'{original_ss_dir}/sources.json'):
        with open(f'{original_ss_dir}/sources.json', encoding='utf-8') as json_file:
            mysrc_nwm_fid = json.load(json_file)
        for key, value in mysrc_nwm_fid.items():
            mysrc_nwm_fid[key] = int(value[0])
    else:
        raise FileNotFoundError('no ele_fid_dict or sources.json found')

    source_eles = Prop(f'{original_ss_dir}/source_sink.in', 2).ip_group[0]
    my_th = TimeHistory(f'{original_ss_dir}/vsource.th', start_time_str=start_time_str, mask_val=None)
    end_time_str = my_th.df.iloc[-1, 0].strftime("%m/%d/%Y, %H:%M:%S")

    my_sources = np.empty((len(source_eles)), dtype=Vsource)
    for i, source_ele in enumerate(source_eles):
        my_sources[i] = Vsource(
            xyz=[hg.xctr[source_ele-1], hg.yctr[source_ele-1], hg.dpe[source_ele-1]],
            hgrid_ie=source_ele,
            nwm_fid=mysrc_nwm_fid[str(source_ele)],
            df=pd.DataFrame({'datetime': my_th.datetime, 'Data': my_th.data[:, i]}),
        )
    print(f'---------------loading vsource took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # ------search USGS stations for each vsource------
    # -----------------------------------------------------------------------------------
    t1 = time.time()
    for i, my_source in enumerate(my_sources):
        # recursively go through upstream branches to find available USGS stations
        find_usgs_along_nwm(
            iupstream=True, starting_seg_id=my_source.nwm_fid,
            vsource=my_source, total_seg_length=0.0, order=0, nwm_shp=nwm_shp
        )
        my_source.upstream_path.writer(
            f'{output_dir}/Diag/Ele{str(my_source.hgrid_ie)}_upstream', shp_gpd=nwm_shp)

        # disable downstream search for now, because it may double count a USGS flow
        # when multiple upstream branches merge into one downstream branch
        # recursively go through downstream branches to find available USGS stations
        # find_usgs_along_nwm(
        #     iupstream=False, starting_seg_id=my_source.nwm_fid,
        #     vsource=my_source, total_seg_length=0.0, order=0, nwm_shp=nwm_shp
        # )
        # my_source.downstream_path.writer(f'{output_dir}/Ele{str(my_source.hgrid_ie)}_downstream', shp_gpd=nwm_shp)

    usgs_stations = [st.st_id for my_source in my_sources for st in my_source.usgs_st]
    usgs_stations = sorted(list(set(usgs_stations)))  # remove duplicates
    print(f'---------------searching usgs stations took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    # Download USGS data
    # -----------------------------------------------------------------------------------
    t1 = time.time()

    # pad one day before and after the start and end time
    padded_start_time = pd.to_datetime(start_time_str) - pd.Timedelta('1 day')
    padded_end_time = pd.to_datetime(end_time_str) + pd.Timedelta('1 day')
    param_dict = {
        "param_id": usgs_var_dict['streamflow']['id'],
        "station_ids": sorted(usgs_stations),
        "cache_fname": (
            f'{cache_folder}/usgs_stofs3d_atl_{usgs_var_dict["streamflow"]["id"]}_'
            f'{padded_start_time.strftime("%Y%m%d")}_{padded_end_time.strftime("%Y%m%d")}.pq'
        ),
        "padded_start_time": padded_start_time.strftime("%Y-%m-%d %H:%M:%S"),
        "padded_end_time": padded_end_time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    param_save = f'{output_dir}/Diag/usgs_download_params1.json'
    with open(param_save, 'w', encoding='utf-8') as f:
        json.dump(param_dict, f, indent=4)
    usgs_data = download_stations(
        param_id=param_dict["param_id"],
        station_ids=param_dict["station_ids"],
        cache_fname=param_dict["cache_fname"],
        datelist=pd.date_range(start=param_dict["padded_start_time"], end=param_dict["padded_end_time"]),
    )
    for data in usgs_data:
        # check data gap
        gap = detect_data_gap(
            df=data['df'], start=start_time_str, end=end_time_str,
            freq=pd.Timedelta('1 hour'))
        if gap > pd.Timedelta('7 days'):
            print(f'warning: USGS station {data["id"]} has a data gap of {gap}, '
                  f'Removing USGS station {data["id"]} from the list of available stations.\n')
            usgs_data.remove(data)
    my_obs = convert_to_ObsData(usgs_data)

    print(f'---------------collecting usgs data took: {time.time()-t1} s ---------------\n')

    # -----------------------------------------------------------------------------------
    #  read NWM data
    # -----------------------------------------------------------------------------------
    t1 = time.time()
    if os.path.exists("nwm.pq"):
        nwm_df = pd.read_parquet("nwm.pq", engine="pyarrow")
    else:
        # assemble subset of fids
        fid_subset = []
        for my_source in my_sources:
            for st in my_source.usgs_st:
                fid_subset.append(st.nearby_nwm_fid)
        my_date_range = pd.date_range(start_time_str, end_time_str, freq="3H")
        stream_flow = read_nwm_data(
            data_dir=nwm_data_dir,
            fids=fid_subset,
            date_range=my_date_range
        )
        # build df
        nwm_df = pd.DataFrame({'datetime': my_date_range}).join(pd.DataFrame(stream_flow, columns=fid_subset))
        # assign utc time zone
        nwm_df['datetime'] = nwm_df['datetime'].dt.tz_localize('UTC')
        nwm_df = nwm_df.loc[:, ~nwm_df.columns.duplicated()]  # drop duplicate columns
        # nwm_df = nwm_df.set_index('datetime').resample("H").mean()
        # nwm_df = nwm_df.interpolate(method='time')  # fill in missing data

        # set datetime as index
        nwm_df = nwm_df.set_index('datetime')

    print(f'---------------reading nwm data took: {time.time()-t1} s ---------------\n')
    # nwm_df.to_parquet("nwm.pq", engine="pyarrow", compression="zstd")
    # -----------------------------------------------------------------------------------
    # Adjust each vsource
    # -----------------------------------------------------------------------------------
    linestyle1 = ['-b', '-g', '-c', '-m', '-y', '-k']
    linestyle2 = ['-.b', '-.g', '-.c', '-.m', '-.y', '-.k']
    for i, my_source in enumerate(my_sources):

        if len(my_source.usgs_st) == 0:  # no usgs to correct nwm
            print(f'warning: source {i}_{my_source.hgrid_ie} does not have nearby usgs obs')
            continue

        # close figure if there is one and make a new one
        plt.close('all')
        plt.figure()

        # original vsource
        df = my_sources[i].df[start_time_str:end_time_str]  # .resample("D").mean()
        nrec = len(df.index)

        # NWM, USGS, adjustment
        # list usgs stations and nwm segs associated with my_source
        mysrc_usgs_id = []
        mysrc_usgs_idx = []
        mysrc_nwm_fid = []
        for st in my_source.usgs_st:
            if st.st_id in my_obs.fID2id.keys():
                idx = my_obs.fID2id[st.st_id]
                tmp_df = my_obs.stations[idx].df
                if len(tmp_df.index) == 0:
                    print(f'warning: USGS {st.st_id} do not have any data'
                          'in the selected period {start_time_str} {end_time_str}')
                else:
                    # nrec_valid_usgs = len(tmp_df.resample("H").mean())
                    # if nrec_valid_usgs == nrec:  # proceed only if there is enough usgs data
                    if tmp_df.index[0] <= pd.to_datetime(start_time_str, utc=True) and \
                       tmp_df.index[-1] >= pd.to_datetime(end_time_str, utc=True):
                        if st.st_id not in mysrc_usgs_id:  # skip duplicates already processed
                            mysrc_usgs_id.append(st.st_id)
                            mysrc_usgs_idx.append(idx)
                            mysrc_nwm_fid.append(st.nearby_nwm_fid)
                    else:
                        print(f'warning: USGS {st.st_id} do not have enough data'
                              'in the selected period {start_time_str} {end_time_str}')
            else:
                print(f'warning: USGS Station {st.st_id} is not in USGS database')
        if len(mysrc_usgs_id) == 0:  # no usgs to correct nwm
            print(f'warning: source {i}_{my_source.hgrid_ie} does not have nearby usgs obs')
            continue

        this_vsource_nwm_array = np.zeros((nrec, len(mysrc_nwm_fid)), dtype=float)
        # interpolate nwm to vsource time
        for k, nwm_fid in enumerate(mysrc_nwm_fid):
            this_vsource_nwm_array[:, k] = np.interp(
                (df.index - df.index[0]).total_seconds().to_numpy(),
                (nwm_df.index - df.index[0]).total_seconds().to_numpy(),
                nwm_df[nwm_fid].to_numpy().squeeze()
            )

        # interpolate usgs to vsource time,
        # where gaps are filled with nwm (this_vs_nwm_array) values
        # (i.e., no difference between usgs and nwm)
        this_vs_usgs_array = np.zeros((nrec, len(mysrc_usgs_id)), dtype=float)
        for k, (nwm_fid, usgs_idx) in enumerate(zip(mysrc_nwm_fid, mysrc_usgs_idx)):
            usgs_station_df = my_obs.stations[usgs_idx].df  # USGS series (cfs) with DatetimeIndex

            this_vs_usgs_array[:, k] = blend_usgs_to_vsource(
                df.index,
                this_vsource_nwm_array[:, k],
                usgs_station_df,
                value_col=1,          # second column is discharge
            )

        plt.figure(figsize=(10, 6))
        plt.plot(df.index, df['Data'], 'k', label="original vsource")
        for k, usgs_id in enumerate(mysrc_usgs_id):
            plt.plot(
                df.index, this_vs_usgs_array[:, k],
                linestyle1[k], linewidth=0.5, label=f'usgs {usgs_id}'
            )
            plt.plot(
                df.index, this_vsource_nwm_array[:, k],
                linestyle2[k], linewidth=0.5, label=f'NWM near usgs {usgs_id}'
            )

        # take the diff (including all relavant usgs-nwm pairs found)
        diff_usgs_nwm = this_vs_usgs_array - this_vsource_nwm_array
        # use maximum to prevent negative sources
        adjusted_source = np.maximum(df['Data'].to_numpy() + np.sum(diff_usgs_nwm, axis=1), 0.0)
        plt.plot(df.index, adjusted_source, '--r', linewidth=1, label="vsource adjusted by USGS-NWM-diff")

        plt.legend()
        # plt.show()
        plt.savefig(f'{output_dir}/Svgs/{i}_{my_source.hgrid_ie}.svg')
        plt.clf()

        # replace vsource at this source
        my_th.df.iloc[:, i+1] = adjusted_source

        print(f'-----------successfully adjusted vsource {i}_{my_source.hgrid_ie}-------------\n')

    my_th.df_propagate()
    my_th.writer(f'{output_dir}/adjusted_vsource.th')


def sample_LA_domain():
    """Sample LA domain for testing purposes."""

    # sample usage (the final product is "adjusted_vsource.th" in the output_dir):

    source_nwm2usgs(
        start_time_str="2024-03-05 00:00:00",
        states=['LA', 'MS', 'TX'],
        f_shapefile="/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2_RiverMapper/shapefiles/LA_nwm_v1p2.shp",
        original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I19/Source_sink/original_source_sink/',
        nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/I09/Source_sink/original_source_sink/20240305/',
        output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I19/Source_sink/USGS_adjusted_sources/',
    )

    # source_nwm2usgs(
    #     start_time_str="2021-08-01 00:00:00",
    #     f_shapefile="/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2_RiverMapper/shapefiles/LA_nwm_v1p2.shp",
    #     original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I09g/Source_sink/original_source_sink/',
    #     nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/I09g/Source_sink/original_source_sink/20210801/',
    #     output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I09j/Source_sink/USGS_adjusted_sources/',
    # )

    # source_nwm2usgs(
    #     start_time_str="2017-12-01 00:00:00",
    #     f_shapefile="/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2_RiverMapper/shapefiles/LA_nwm_v1p2.shp",
    #     original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I14/Source_sink/USGS_adjusted_sources/original/',
    #     nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/I13/Source_sink/original_source_sink/20171201/',
    #     output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I14/Source_sink/USGS_adjusted_sources/',
    # )

    # source_nwm2usgs(
    #     start_time_str="2017-12-01 00:00:00",
    #     f_shapefile="/sciclone/schism10/Hgrid_projects/NWM/ecgc/ecgc.shp",
    #     original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15_v7/Source_sink/original_source_sink/',
    #     nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/I13/Source_sink/original_source_sink/20171201/',
    #     output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15_v7/Source_sink/USGS_adjusted_sources/',
    # )


def sample_stofs_v7v8():
    """
    Sample STOFS v7/v8 domain for testing purposes.
    """
    source_nwm2usgs(
        start_time_str="2017-12-01 00:00:00",
        states=STOFS3D_ATL_STATES,
        f_shapefile="/sciclone/schism10/Hgrid_projects/NWM/ecgc/ecgc.shp",
        original_ss_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15n_v7/original_source_sink0/',
        nwm_data_dir='/sciclone/schism10/feiye/STOFS3D-v8/NWM/CONUS/netcdf/CHRTOUT/for_2018_hindcast/',
        output_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15n_v7/Source_sink/USGS_adjusted_sources/',
    )


if __name__ == "__main__":
    sample_LA_domain()
    print('done')
