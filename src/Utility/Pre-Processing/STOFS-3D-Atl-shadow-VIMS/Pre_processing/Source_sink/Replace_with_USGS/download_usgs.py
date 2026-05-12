# -*- coding: utf-8 -*-
"""
Requires:
    pip install dataretrieval pandas numpy gsw requests

Notes:
- For "state discovery", NWIS site service may not always include begin/end dates in the
  response; we keep columns for compatibility but leave them empty if not provided.
- IV service returns UTC timestamps. If you need local time, convert explicitly later.
"""

import os
import pickle
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import time
from requests.exceptions import SSLError, ReadTimeout, ConnectionError as ReqConnError

import gsw

from dataretrieval import nwis
# from schism_py_pre_post.Download.Data import ObsData, Station   # keep your original import
# If needed, uncomment the above; left commented here to make this module self-contained.


# -----------------------------
# Configuration & helpers
# -----------------------------

usgs_var_dict = {
    "streamflow":  {"id": "00060", "unit_conv": 0.028316846592},  # cfs -> m^3/s
    "salinity":    {"id": "00480", "unit_conv": 1},               # PSU directly (rarely in NWIS)
    "gauge height":{"id": "00065", "unit_conv": 1},               # ft
    "temperature": {"id": "00010", "unit_conv": 1},               # degC
    "conductance": {"id": "00095", "unit_conv": 0.001},           # uS/cm -> mS/cm
}

ecgc_states = [
    'ME','NH','MA','RI','CT','NY','NJ','DE','PA','MD','DC','VA','NC','SC',
    'GA','FL','AL','MI','LA','TX'
]


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


# -----------------------------
# Data containers (unchanged)
# -----------------------------

class GenericObsData:
    __slots__ = ['station_info', 'df']
    def __init__(self, station_info=None, df=None):
        self.station_info = station_info
        self.df = df


# -----------------------------
# Station discovery (states)
# -----------------------------

def get_usgs_stations_from_state(states=('VA', 'MD'), data_type='iv', parameter=None):
    """
    Discover USGS stations by states using NWIS 'site' service via dataretrieval.

    Parameters
    ----------
    states : iterable of state abbreviations
    data_type : 'iv' or 'dv' (filters on hasDataTypeCd)
    parameter : optional parameter code to filter, e.g. '00060'

    Returns
    -------
    DataFrame with columns:
      ['site_no','station_nm','state','dec_long_va','dec_lat_va','begin_date','end_date']
    begin_date/end_date may be NaT if not provided by the API response.
    """
    assert data_type in ('iv','dv'), "data_type must be 'iv' or 'dv'"
    rows = []
    for state in states:
        # Build filters. dataretrieval passes **kwargs to NWIS site API.
        # https://waterservices.usgs.gov/rest/Site-Service.html
        kw = dict(
            stateCd=state,
            siteStatus='all',
            hasDataTypeCd=data_type
        )
        if parameter:
            kw['parameterCd'] = parameter

        df_info, _meta = nwis.get_info(**kw)

        if df_info is None or len(df_info) == 0:
            print(f"{state}: 0 stations found.")
            continue

        # Standardize columns we care about (guard against missing fields)
        keep = {
            'site_no': 'site_no',
            'station_nm': 'station_nm',
            'state_cd': 'state',           # rename to 'state'
            'dec_long_va': 'dec_long_va',
            'dec_lat_va': 'dec_lat_va',
            # begin_date/end_date may exist depending on filters; handle gracefully
            'begin_date': 'begin_date',
            'end_date': 'end_date',
            # Sometimes 'parm_cd' column appears; we don't strictly need it for output
        }
        df_keep = pd.DataFrame()
        for src, dst in keep.items():
            if src in df_info.columns:
                df_keep[dst] = df_info[src]
            else:
                # ensure column exists
                df_keep[dst] = pd.NaT if 'date' in src else None

        print(f"{state} stations found: {len(df_keep)}")
        rows.append(df_keep)

    if rows:
        out = pd.concat(rows, ignore_index=True)
    else:
        out = pd.DataFrame(columns=['site_no','station_nm','state','dec_long_va','dec_lat_va','begin_date','end_date'])

    print(f"Total stations found: {len(out)}")
    return out


# -----------------------------
# Core download via dataretrieval
# -----------------------------

def _discover_station_metadata(station_ids):
    """
    Pull station name/lat/lon for a list of sites with one NWIS 'site' query.
    Returns a dict: site_no -> dict(name, lon, lat)
    """
    if isinstance(station_ids, (list, tuple)):
        sites = ",".join(station_ids)
    else:
        sites = str(station_ids)
    df_info, _ = nwis.get_info(sites=sites, siteStatus='all')
    name_map = {}
    if df_info is not None and len(df_info) > 0:
        for _, r in df_info.iterrows():
            name_map[str(r['site_no'])] = {
                'name': r.get('station_nm', ''),
                'lon': float(r.get('dec_long_va', np.nan)) if not pd.isna(r.get('dec_long_va', np.nan)) else np.nan,
                'lat': float(r.get('dec_lat_va', np.nan)) if not pd.isna(r.get('dec_lat_va', np.nan)) else np.nan,
            }
    return name_map


def _extract_value_col(df, parameter):
    """
    Given a dataretrieval iv/dv dataframe reset to columns,
    find the value column for parameter code, plus qualifier column if present.
    """
    p = str(parameter)
    cand = [c for c in df.columns if c.startswith(p)]
    # Prefer explicit suffixes
    value_col = None
    if f"{p}_Value" in cand:
        value_col = f"{p}_Value"
    elif f"{p}_Mean" in cand:
        value_col = f"{p}_Mean"
    elif p in df.columns:
        value_col = p
    else:
        # fall back: pick first that starts with p and does not end with '_cd'
        tmp = [c for c in cand if not c.endswith('_cd')]
        if tmp:
            value_col = tmp[0]
    qual_col = f"{p}_cd" if f"{p}_cd" in df.columns else None
    return value_col, qual_col


def _site_has_param(site_id, param_id, data_type='iv'):
    # Uses site service with series catalog; fast and avoids blind retries
    try:
        df, _ = nwis.get_info(sites=str(site_id), siteStatus='all',
                              seriesCatalogOutput=True)
        if df is None or 'series_catalog' not in df.columns:
            return True  # don't block if catalog missing
        sc = df['series_catalog'].iloc[0]
        if isinstance(sc, list):
            return any(str(param_id) in str(item) and data_type in str(item).lower() for item in sc)
        return True
    except Exception:
        return True


def _get_iv_retry(sites_arg, param_id, start, end, tries=3, backoff=1.5):
    for k in range(tries):
        try:
            return nwis.get_iv(sites=sites_arg, parameterCd=param_id,
                               start=start, end=end, siteStatus='all')
        except (SSLError, ReadTimeout, ReqConnError) as e:
            wait = backoff * (2**k)
            print(f"get_iv network error: {e} → retry in {wait:.1f}s")
            time.sleep(wait)
        except Exception as e:
            print(f"get_iv non-retryable error: {e}")
            return (None, None)
    return (None, None)


def download_stations(
    param_id=None, station_ids=None, cache_fname=None,
    datelist=pd.date_range(start='2000-01-01', end='2000-01-02'),
    stations_chunk_size=50
):
    """
    Download instantaneous values (IV) for a list of stations and a date range
    using dataretrieval.nwis.get_iv, returning a list of GenericObsData.

    (Matches your original download_stations() return format.)
    """
    if station_ids is None:
        station_ids = []
    if isinstance(station_ids, np.ndarray):
        station_ids = station_ids.tolist()

    if cache_fname is not None and os.path.exists(cache_fname):
        print(f'Loading cached data from {cache_fname} ...')
        with open(cache_fname, 'rb') as f:
            total_data = pickle.load(f)
        return total_data

    start = pd.to_datetime(datelist[0]).strftime('%Y-%m-%d')
    end   = pd.to_datetime(datelist[-1]).strftime('%Y-%m-%d')

    print(f'Downloading data for {len(station_ids)} stations ...')
    print(f'station ids: {station_ids}\n')

    # Pre-fetch station metadata (name/lat/lon)
    meta_map = _discover_station_metadata(station_ids)

    total_data = []
    for i, station_ids_chunk in enumerate(chunks(station_ids, stations_chunk_size)):
        sites_chunk = [str(s) for s in station_ids_chunk]

        # 1) multi-site attempt with retries
        df_iv, _ = _get_iv_retry(",".join(sites_chunk), str(param_id), start, end, tries=4, backoff=1.5)
        if df_iv is None:
            df_iv = pd.DataFrame()

        if len(df_iv) > 0:
            df_iv = df_iv.reset_index()

        value_col, _qual = (None, None)
        if len(df_iv) > 0:
            value_col, _qual = _extract_value_col(df_iv, str(param_id))

        # Figure out which sites are missing
        got_sites = set(df_iv['site_no'].astype(str)) if len(df_iv) > 0 else set()
        missing_sites = [s for s in sites_chunk if s not in got_sites]

        # 2) per-site fallback for missing ones (with retries)
        fallback_frames = []
        for s in missing_sites:
            if not _site_has_param(s, param_id, 'iv'):
                # Skip silent no-data stations quickly
                print(f"  Skip {s}: no IV series for {param_id} in catalog")
                continue
            df1, _ = _get_iv_retry(s, str(param_id), start, end, tries=4, backoff=1.5)
            if df1 is not None and len(df1) > 0:
                fallback_frames.append(df1.reset_index())

        if fallback_frames:
            extra = pd.concat(fallback_frames, ignore_index=True, sort=False)
            if len(df_iv) == 0:
                df_iv = extra
            else:
                # align columns
                common = list(set(df_iv.columns).intersection(set(extra.columns)))
                df_iv = pd.concat([df_iv[common], extra[common]], ignore_index=True)

        if df_iv is None or len(df_iv) == 0:
            print(f"Chunk {i}: no IV data returned.")
            print(f'stations processed: {min(len(station_ids),(i+1)*stations_chunk_size)} of {len(station_ids)}')
            continue

        # recompute value_col after concat if needed
        if value_col is None or value_col not in df_iv.columns:
            value_col, _qual = _extract_value_col(df_iv, str(param_id))
            if value_col is None:
                print(f"Chunk {i}: cannot locate value column for parameter {param_id}. "
                      f"Columns: {list(df_iv.columns)}")
                print(f'stations processed: {min(len(station_ids),(i+1)*stations_chunk_size)} of {len(station_ids)}')
                continue

        # Build GenericObsData per site
        for site, g in df_iv.groupby('site_no', as_index=False):
            # g has at least ['datetime', value_col]
            dates = pd.to_datetime(g['datetime'], utc=True)
            values = pd.to_numeric(g[value_col], errors='coerce')

            df = pd.DataFrame({'date': dates, 'value': values})

            # station meta
            meta = meta_map.get(str(site), {'name':'', 'lon':np.nan, 'lat':np.nan})
            station_info = {
                'id': str(site),
                'name': meta.get('name', ''),
                'lon': meta.get('lon', np.nan),
                'lat': meta.get('lat', np.nan),
                'var_name': param_id,      # keep same fields you expected before
                'var_code': param_id,
                'unit': ''                 # unit not always exposed here; can be added via pmcodes if needed
            }

            total_data.append(GenericObsData(station_info=station_info, df=df))

        print(f'stations processed: {min(len(station_ids),(i+1)*stations_chunk_size)} of {len(station_ids)}')

    print(f'Number of stations with available data: {len(total_data)}')

    # save to cache
    if cache_fname is not None:
        print(f'Saving downloaded data to {cache_fname} ...')
        with open(cache_fname, 'wb') as f:
            pickle.dump(total_data, f)

    return total_data


# -----------------------------
# Utilities: averaging & writers
# -----------------------------

def write_time_average(input_data, param_id=None, unit_conv=1, outfilename=None):
    """
    Write mean values as 'x y z id' lines.
    If param is '00095' (conductance), convert to salinity (PSU) at 25C, 0 dbar.
    """
    if outfilename is None:
        raise ValueError("outfilename must be provided")

    with open(outfilename, 'w') as fout:
        for data in input_data:
            x = data.station_info['lon']
            y = data.station_info['lat']
            z = data.df['value'].mean()
            if param_id == '00095':  # convert conductance (mS/cm after unit_conv) to salinity
                z = gsw.conversions.SP_from_C(z * unit_conv, 25, 0)
            if np.isnan(z):
                continue
            fout.write(f"{x} {y} {z} {data.station_info['id']}\n")


def all_states_time_average(param_id=None, unit_conv=1, states=None,
                            datelist=pd.date_range(start='2000-01-01', end='2000-01-02'),
                            outfilename=None):
    """
    Discover stations by states and write mean XYZ for the date window.
    Uses NWIS site service + iv pulls (dataretrieval).
    """
    if outfilename is None:
        raise ValueError("outfilename must be provided")

    if os.path.exists(outfilename):
        print(f'{outfilename} exists, skipping ...')
        return

    if states is None:
        states = ecgc_states

    # Discover stations for the requested parameter (faster/lower-noise)
    station_info_df = get_usgs_stations_from_state(states, data_type='iv', parameter=param_id)
    station_ids = station_info_df["site_no"].astype(str).tolist()

    Path(os.path.dirname(outfilename) or ".").mkdir(parents=True, exist_ok=True)

    total_data = download_stations(
        param_id=param_id,
        station_ids=station_ids,
        datelist=datelist
    )
    write_time_average(
        input_data=total_data,
        param_id=param_id,
        unit_conv=unit_conv,
        outfilename=outfilename
    )


# -----------------------------
# Single-station helper
# -----------------------------

def download_single_station(
    station_id='04044755',
    param_id='00010',
    datelist=pd.date_range(start='2019-10-10', end='2019-11-10')
):
    """
    Minimal clone of your single-station IV pull returning [GenericObsData]
    """
    start = pd.to_datetime(datelist[0]).strftime('%Y-%m-%d')
    end   = pd.to_datetime(datelist[-1]).strftime('%Y-%m-%d')

    df_iv, _ = nwis.get_iv(
        sites=str(station_id),
        parameterCd=str(param_id),
        start=start,
        end=end,
        siteStatus='all'
    )
    if df_iv is None or len(df_iv) == 0:
        return []

    df_iv = df_iv.reset_index()
    value_col, _qual = _extract_value_col(df_iv, str(param_id))
    if value_col is None:
        return []

    dates = pd.to_datetime(df_iv['datetime'], utc=True)
    values = pd.to_numeric(df_iv[value_col], errors='coerce')
    df = pd.DataFrame({'date': dates, 'value': values})

    # station meta
    meta_map = _discover_station_metadata([station_id])
    meta = meta_map.get(str(station_id), {'name':'', 'lon':np.nan, 'lat':np.nan})
    station_info = {
        'id': str(station_id),
        'name': meta.get('name', ''),
        'lon': meta.get('lon', np.nan),
        'lat': meta.get('lat', np.nan),
        'var_name': param_id,
        'var_code': param_id,
        'unit': ''
    }
    return [GenericObsData(station_info=station_info, df=df)]


# -----------------------------
# (Optional) Old ObsData format
# -----------------------------

def convert_to_ObsData(total_data, cache_fname=None):
    """
    Convert the downloaded "total_data" (list of GenericObsData) to your older ObsData format.
    """
    from schism_py_pre_post.Download.Data import ObsData, Station  # local import to avoid hard dependency here

    stations = []
    for data in total_data:
        my_station = Station(fname=None)
        my_station.id = data.station_info['id']
        my_station.lon = data.station_info['lon']
        my_station.lat = data.station_info['lat']
        my_station.df = data.df.set_index(pd.DatetimeIndex(pd.to_datetime(data.df["date"], utc=True)))
        stations.append(my_station)

    obs_data = ObsData(fname=None)
    obs_data.stations = stations
    obs_data.set_fID2id()

    if cache_fname is not None:
        obs_data.saved_file = cache_fname
        obs_data.save(obs_data.saved_file)

    return obs_data


# -----------------------------
# Example: station info sample
# -----------------------------

def sample_get_station_info():
    """
    Example usage of get_usgs_stations_from_state + a targeted fetch list,
    mirroring your original function shape.
    """
    stations = [
        "02492511", "02492519", "02492600", "02492700", "07374000", "073745245", "073745253", "07374581",
        "07375050", "07375170", "07375175", "07375222", "07375230", "07375500", "07375650", "07376000",
        "07378050", "07378500", "07378745", "07378746", "07378748", "07378810", "07379050", "07379075",
        "07380101", "07380102", "07380120", "07380126", "07380127", "07380200", "07380212", "07380215",
        "073802220", "073802225", "0738022295", "0738022395", "073802245", "073802273", "073802280", "073802282",
        "073802284", "07380330", "07380401", "07380500", "07381000", "07381150", "07381324", "07381350",
        "07381355", "07381450", "07381454", "07381460", "07381490", "07381515", "073815450", "07381590",
        "073815945", "073815963", "07381600", "07384400", "07385700", "07385702", "07385765", "07385820",
        "07386600", "07386850", "08010000", "08012150", "292952090565300", "293809092361500", "294045092492300",
        "294717092250000", "295011091184300", "2951190901217", "295124089542100", "295447091191500",
        "295501090190400", "300312091320000", "300507091355600", "300602090375100",
        "300703089522700", "301200090072400", "301324090382400", "302020091435700",
    ]

    total_data = download_stations(
        param_id=usgs_var_dict['gauge height']['id'],
        station_ids=stations,
        datelist=pd.date_range(start='2021-07-31', end='2021-08-01')
    )

    # Save station id and names in CSV (dedup by site_no)
    station_ids = [data.station_info['id'] for data in total_data]
    station_names = [data.station_info['name'] for data in total_data]
    station_ids = [sid.strip().strip('"') for sid in station_ids]
    station_names = [sname.strip().strip('"') for sname in station_names]

    df = pd.DataFrame({'site_no': station_ids, 'station_nm': station_names}).drop_duplicates(subset='site_no')
    df.to_csv('station_info.csv', index=False, sep=';')
    print(df.head())


# -----------------------------
# Main
# -----------------------------

if __name__ == "__main__":
    sample_get_station_info()
