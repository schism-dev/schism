# -*- coding: utf-8 -*-
# USGS IV downloader without climata (RDB-based)

from __future__ import annotations
import os
import re
import json
from typing import List, Iterable, Tuple, Optional
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

import gsw  # for conductance->salinity conversion if desired

# If you rely on these types elsewhere:
from schism_py_pre_post.Download.Data import ObsData, Station


# -------------------- Configuration --------------------

# dict of param ids and unit conversions (to your preferred internal units)
usgs_var_dict = {
    "streamflow":   {"id": "00060", "unit_conv": 0.028316846592},  # cfs -> m^3/s
    "salinity":     {"id": "00480", "unit_conv": 1.0},             # PSU (rarely available)
    "gauge height": {"id": "00065", "unit_conv": 1.0},             # feet or ft? (check site metadata)
    "temperature":  {"id": "00010", "unit_conv": 1.0},             # degC
    "conductance":  {"id": "00095", "unit_conv": 0.001},           # uS/cm -> mS/cm (for SP_from_C)
}

# East & Gulf Coast states default (you can override)
ecgc_states = ['ME', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'DE', 'PA', 'MD', 'DC',
               'VA', 'NC', 'SC', 'GA', 'FL', 'AL', 'MI', 'LA', 'TX']

USGS_IV_URL = "https://waterservices.usgs.gov/nwis/iv"
USGS_SITE_URL = "https://waterservices.usgs.gov/nwis/site"


# ------------------------Caching------------------------
class GenericObsData:
    __slots__ = ['station_info', 'df']

    def __init__(self, station_info=None, df=None):
        self.station_info = station_info or {}
        self.df = df  # expects columns: ['date','value', (optional) 'tz_cd']


# ----- Flatten to a single DataFrame -----
def to_portable_df(objs: List[GenericObsData]) -> pd.DataFrame:
    rows = []
    for o in objs:
        info = o.station_info
        df = o.df.copy()

        # Ensure a 'date' column (if it's the index, reset it)
        if df.index.name == "date":
            df = df.reset_index()

        # Normalize datetime to UTC (portable and unambiguous)
        if "date" in df.columns:
            if isinstance(df['date'].dtype, pd.DatetimeTZDtype):
                df["date"] = df["date"].dt.tz_convert("UTC")
            else:
                raise ValueError("DataFrame 'date' column must be timezone-aware datetime")

        site_no = str(info.get("id", ""))
        station_nm = info.get("name", None)
        lon = info.get("lon", None)
        lat = info.get("lat", None)
        var_code = info.get("var_code", None)

        block = df.loc[:, ["date", "value"]].copy()
        block.insert(0, "site_no", site_no)
        block.insert(1, "station_nm", station_nm)
        block.insert(2, "lon", lon)
        block.insert(3, "lat", lat)
        block.insert(4, "var_code", var_code)
        rows.append(block)

    if not rows:
        return pd.DataFrame(columns=["site_no", "station_nm", "lon", "lat", "var_code", "date", "value"])

    out = pd.concat(rows, ignore_index=True)

    # Enforce types
    out["site_no"] = out["site_no"].astype(str)
    out["date"] = pd.to_datetime(out["date"], utc=True, errors="coerce")
    out["lon"] = pd.to_numeric(out["lon"], errors="coerce")
    out["lat"] = pd.to_numeric(out["lat"], errors="coerce")
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    return out


# ----- Save as Parquet (portable) -----
def save_cache_parquet(objs: List[GenericObsData], path: str) -> None:
    df = to_portable_df(objs)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, engine="pyarrow", compression="zstd", index=False)


# Optional CSV fallback (slower, bigger, but universal)
def save_cache_csv(objs: List[GenericObsData], path: str) -> None:
    df = to_portable_df(objs)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


# ----- Read and (optionally) reconstruct your objects -----
def load_cache_parquet(path: str) -> pd.DataFrame:
    return pd.read_parquet(path, engine="pyarrow")


def df_to_generic_obs(df: pd.DataFrame) -> List[GenericObsData]:
    out: List[GenericObsData] = []
    for site_no, g in df.groupby("site_no", sort=False):
        info = {
            "id": site_no,
            "name": g["station_nm"].iloc[0] if "station_nm" in g.columns else None,
            "lon": float(g["lon"].iloc[0]) if "lon" in g.columns else None,
            "lat": float(g["lat"].iloc[0]) if "lat" in g.columns else None,
            "var_name": g["var_code"].iloc[0] if "var_code" in g.columns else None,
            "var_code": g["var_code"].iloc[0] if "var_code" in g.columns else None,
            "unit": None,
        }
        sub = g[["date", "value"]].copy()
        # keep tz-aware UTC datetimes as-is
        sub = sub.reset_index(drop=True)
        out.append(GenericObsData(station_info=info, df=sub))
    return out


def _normalize_stations(station_ids):
    # keep digits only / or just strip if you already cleaned earlier; here we just coerce to str and strip
    ids = [str(s).strip() for s in station_ids if s is not None and str(s).strip() != ""]
    # de-duplicate while preserving original relative order, then sort for canonical form
    seen = set()
    uniq = []
    for s in ids:
        if s not in seen:
            uniq.append(s)
            seen.add(s)
    return sorted(uniq)  # canonical order


def write_manifest(path: str, param_id: str, start, end, station_ids):
    stations = _normalize_stations(station_ids)
    start_iso = pd.to_datetime(start, utc=True).isoformat().replace("+00:00", "Z")
    end_iso = pd.to_datetime(end, utc=True).isoformat().replace("+00:00", "Z")

    manifest = {
        "version": 1,
        "parameter": param_id,
        "start": start_iso,
        "end": end_iso,
        "stations": sorted(stations),       # explicit list (nice for debugging)
        "count": len(stations),
    }
    mpath = Path(path).with_suffix(Path(path).suffix + ".json")  # e.g., data.parquet.json
    mpath.parent.mkdir(parents=True, exist_ok=True)
    with open(mpath, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True, ensure_ascii=False)


def check_manifest(
    cache_path: str, station_ids: list[str],
    param_id: str | None = None, start=None, end=None
) -> bool:

    mpath = Path(cache_path).with_suffix(Path(cache_path).suffix + ".json")
    if not mpath.exists():
        print("Manifest missing.")
        return False

    with open(mpath, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    # version
    if manifest.get("version") != 1:
        print("Manifest version mismatch!")
        return False

    # stations
    want = _normalize_stations(station_ids)
    have = manifest.get("stations", [])
    if want != have:
        print("Station list mismatch!")
        return False

    # parameter (optional check if provided)
    if param_id is not None and manifest.get("parameter") != param_id:
        print("Parameter mismatch!")
        return False

    # time window (optional checks if provided)
    def _to_iso_z(x):
        return pd.to_datetime(x, utc=True).isoformat().replace("+00:00", "Z")
    if start is not None and manifest.get("start") != _to_iso_z(start):
        print("Start time mismatch!")
        return False
    if end is not None and manifest.get("end") != _to_iso_z(end):
        print("End time mismatch!")
        return False

    print("Manifest OK.")
    return True


# -------------------- Session / Retry --------------------

def make_session(ignore_proxies: bool = True) -> requests.Session:
    """Requests session with retries, backoff, and good pool sizes."""
    retry = Retry(
        total=6, connect=3, read=3, status=3,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods={"GET"},
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=32, pool_maxsize=32)
    s = requests.Session()
    s.mount("https://", adapter)
    s.headers.update({
        "User-Agent": "VIMS-STOFS usgs-iv-downloader (requests)",
        "Accept": "text/plain",
    })
    # ignore system proxy env vars unless you truly need them
    s.trust_env = not ignore_proxies
    return s


# -------------------- Utilities --------------------

def chunks(lst: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def chunk_range(start: datetime, end: datetime, days: int) -> Iterable[Tuple[datetime, datetime]]:
    cur = start
    step = timedelta(days=days)
    while cur < end:
        nxt = min(cur + step, end)
        yield cur, nxt
        cur = nxt


# -------------------- RDB Parsing --------------------

def parse_rdb(text: str) -> pd.DataFrame:
    """
    Parse USGS RDB into DataFrame.
    Skips comment lines starting with '#'.
    The first non-comment line is header; second non-comment line is dtypes; then data rows.
    Returns a DataFrame with string columns; cast as needed downstream.
    """
    lines = text.splitlines()
    # Strip comment lines
    non_comment = [ln for ln in lines if ln and not ln.startswith("#")]

    if len(non_comment) < 3:
        # header, dtype, maybe no data
        return pd.DataFrame()

    header = non_comment[0].split('\t')
    # dtype line:
    # dtype = non_comment[1].split('\t')  # not strictly needed usually
    data_lines = non_comment[2:]

    rows = []
    for ln in data_lines:
        if not ln.strip():
            continue
        parts = ln.split('\t')
        # Ensure length match (USGS sometimes leaves trailing tabs)
        if len(parts) < len(header):
            parts += [""] * (len(header) - len(parts))
        rows.append(parts[:len(header)])

    df = pd.DataFrame(rows, columns=header)
    return df


# -------------------- Site metadata --------------------

def usgs_site_info(session: requests.Session, site_ids: List[str]) -> pd.DataFrame:
    """
    Fetch site metadata (name, lat/lon, etc.) for a list of sites via RDB.
    Returns DataFrame at least with columns: site_no, station_nm, dec_long_va, dec_lat_va
    """
    # USGS allows many sites per request, but stay conservative (e.g., 200/site batch)
    out = []
    for group in chunks(site_ids, 200):
        params = {
            "format": "rdb",
            "sites": ",".join(group),
            "siteStatus": "all",
            "siteOutput": "Expanded",
        }
        r = session.get(USGS_SITE_URL, params=params, timeout=(5, 120))
        r.raise_for_status()
        df = parse_rdb(r.text)
        if not df.empty:
            out.append(df)

    if not out:
        return pd.DataFrame(columns=["site_no", "station_nm", "dec_long_va", "dec_lat_va"])

    meta = pd.concat(out, ignore_index=True).drop_duplicates(subset=["site_no"])
    # Ensure numeric types for coords if present
    for c in ("dec_long_va", "dec_lat_va"):
        if c in meta.columns:
            meta[c] = pd.to_numeric(meta[c], errors="coerce")
    return meta


def get_usgs_stations_from_state(states: List[str] = None,
                                 data_type: str = 'iv',
                                 parameter: Optional[str] = None) -> pd.DataFrame:
    """
    List stations by state(s) using the site service (RDB).
    Returns columns: site_no, station_nm, state, dec_long_va, dec_lat_va, begin_date, end_date

    Note:
      - We filter by hasDataTypeCd=iv when data_type='iv'
      - If parameter is provided, also pass parameterCd filter (narrower & faster)
    """
    if states is None:
        states = ecgc_states

    session = make_session(ignore_proxies=True)
    collected = []

    for st in states:
        params = {
            "format": "rdb",
            "stateCd": st,
            "siteStatus": "all",
            "siteOutput": "Expanded",
        }
        if data_type == "iv":
            params["hasDataTypeCd"] = "iv"
        if parameter:
            params["parameterCd"] = parameter

        # Retry a few times here manually (we already have Retry at session level)
        for attempt in range(5):
            try:
                r = session.get(USGS_SITE_URL, params=params, timeout=(5, 120))
                r.raise_for_status()
                df = parse_rdb(r.text)
                if not df.empty:
                    df["state"] = st
                    collected.append(df)
                print(f"{st} stations found: {0 if df is None else len(df)}")
                break
            except Exception as e:
                print(f"State {st}: attempt {attempt+1}/5 failed: {e}")
                if attempt == 4:
                    raise

    if not collected:
        return pd.DataFrame(columns=[
            'site_no', 'station_nm', 'state', 'dec_long_va', 'dec_lat_va', 'begin_date', 'end_date'])

    big = (pd.concat(collected, ignore_index=True)
             .drop_duplicates(subset=["site_no"]))

    # Keep the columns you used before (create if missing)
    cols = ['site_no', 'station_nm', 'state', 'dec_long_va', 'dec_lat_va', 'begin_date', 'end_date']
    for c in cols:
        if c not in big.columns:  # initialize if not present
            big[c] = np.nan

    # cast coords
    big["dec_long_va"] = pd.to_numeric(big["dec_long_va"], errors="coerce")
    big["dec_lat_va"] = pd.to_numeric(big["dec_lat_va"], errors="coerce")

    print(f"Total stations found: {len(big)}")
    return big[cols]


# -------------------- IV fetch (RDB) --------------------
def fetch_iv_rdb(session: requests.Session,
                 sites: list[str],
                 parameter: str,
                 start: datetime,
                 end: datetime,
                 days_per_chunk: int = 30,
                 data_gap_okay: bool = False,
                 timeout: tuple[float, float] = (5.0, 180.0)) -> pd.DataFrame:
    """
    Fetch IV data for MULTIPLE sites and ONE parameter as RDB.
    Returns tidy DataFrame with columns ['site_no','date','value','tz_cd'].
    Keeps timestamps in local time (no UTC conversion).

    data_gap_okay: If True, allows continuing past failures to a time chunk (only for single sites).
    If multiple sites are requested (len(sites)>1), any failure raises an exception immediately.

    A practical parameter setting is to set days_per_chunk to the tolerable maximum data gap,
    and set data_gap_okay=False, so that a site is skipped upon the first failure to retrieve
    a time chunk, thus avoiding long waits on repeated failures for the same site.

    Example:
        df = fetch_iv_rdb(sess, ["07374000","07374500"], "00060",
                          pd.Timestamp("2020-06-01"), pd.Timestamp("2020-06-05"))
    """
    blocks: list[pd.DataFrame] = []

    if isinstance(sites, str):
        sites = [sites]  # single site to list

    # break long date spans into smaller chunks
    for sdt, edt in chunk_range(start, end, days=days_per_chunk):
        params = {
            "format": "rdb",
            "sites": ",".join(sites),
            "parameterCd": parameter,
            "startDT": sdt.isoformat(timespec="seconds"),
            "endDT": edt.isoformat(timespec="seconds"),
            "siteStatus": "all",
        }

        try:
            r = session.get(USGS_IV_URL, params=params, timeout=timeout)
        except Exception as e:
            print(f"Error fetching data for sites {sites[0]}–{sites[-1]} from {sdt} to {edt}: {e}")
            if len((sites)) > 1:  # only raise if multi-site fetch fails, single-site fetch allows data gaps
                raise e
            else:
                if data_gap_okay:
                    print("Continuing despite error in single site fetch.")
                    continue
                else:
                    print(f"Failed single-site fetch for {sites[0]} during {sdt} to {edt} ({days_per_chunk} days).")
                    print("Giving up on this station.")
                    raise e

        df = parse_rdb(r.text)
        if df.empty:
            continue

        # --- detect all value columns that correspond to the parameter ---
        val_cols = [c for c in df.columns if c.endswith(f"_{parameter}") and not c.endswith("_cd")]
        if parameter in df.columns:
            val_cols.append(parameter)
        if not val_cols:
            continue  # no data in this chunk

        # --- melt each station's column into long format ---
        # IV RDB includes columns: site_no, datetime, tz_cd, {siteid_paramCd}
        # e.g. "07374000_00060", "07374500_00060"
        common_cols = [c for c in ("agency_cd", "site_no", "datetime", "tz_cd") if c in df.columns]
        for col in val_cols:
            sub = df.loc[:, common_cols + [col]].copy()
            sub.rename(columns={col: "value"}, inplace=True)
            # Extract site_no from the column name if possible (some files omit 'site_no' field)
            if "site_no" not in sub.columns or sub["site_no"].isna().all():
                m = re.match(r"(\d{7,})_", col)
                if m:
                    sub["site_no"] = m.group(1)
            sub["value"] = pd.to_numeric(sub["value"], errors="coerce")
            sub = sub.dropna(subset=["datetime", "value"])
            sub["date"] = pd.to_datetime(sub["datetime"], errors="coerce")
            sub = sub.dropna(subset=["date"])
            blocks.append(sub[["site_no", "date", "value", "tz_cd"]])

    if not blocks:
        return pd.DataFrame(columns=["site_no", "date", "value", "tz_cd"])

    out = (
        pd.concat(blocks, ignore_index=True)
        .sort_values(["site_no", "date"])
        .drop_duplicates(subset=["site_no", "date"], keep="last")
        .reset_index(drop=True)
    )
    return out


# -------------------- Main download logic (no climata) --------------------
tz_map = {
    "EST": "US/Eastern", "EDT": "US/Eastern",
    "CST": "US/Central", "CDT": "US/Central",
    "MST": "US/Mountain", "MDT": "US/Mountain",
    "PST": "US/Pacific", "PDT": "US/Pacific",
    "AKST": "US/Alaska", "AKDT": "US/Alaska",
    "HST": "US/Hawaii",
    "UTC": "UTC", "GMT": "UTC",
}


def detect_data_gap(df: pd.DataFrame, start, end, *, freq: pd.Timedelta | None = None) -> pd.Timedelta:
    """
    Return the largest *missing* duration within [start, end] for a time-indexed df.
    - If `freq` is None, infer it from the median spacing of existing timestamps. Note this needs not to be the true dt, but the tolerable gap, e.g., the native dt is 15 min, but if you can accept up to 1 hour gaps, set freq=1H.
    - A gap between two consecutive timestamps t[i-1] -> t[i] contributes
      max(0, (t[i] - t[i-1]) - freq) to the missing time.
    - Also considers leading (start -> first ts) and trailing (last ts -> end) coverage,
      which are counted as-is (not reduced by freq).
    Assumes df.index is a DatetimeIndex (unsorted OK).
    """
    start = pd.to_datetime(start)
    end = pd.to_datetime(end)

    if start >= end:
        return pd.Timedelta(0)

    # Trim to range
    if df.empty:
        return end - start

    sub = df.loc[start:end]
    if sub.empty:
        return end - start

    ts = sub.index.to_series().dropna().sort_values().unique()
    if len(ts) == 0:
        return end - start

    # Infer freq if not supplied
    if freq is None:
        if len(ts) < 2:
            # Can't infer cadence from a single point; whole covered span apart from that one point is unknown.
            # Be conservative: the maximum gap is whichever side is longer.
            return max(ts.iloc[0] - start, end - ts.iloc[-1])
        deltas = pd.Series(ts).diff().dropna()
        # Use median to be robust to outliers/occasional big gaps
        freq = deltas.median()
        # Guard against zero/negative
        if pd.isna(freq) or freq <= pd.Timedelta(0):
            # Fall back to minimum positive delta if available
            pos = deltas[deltas > pd.Timedelta(0)]
            freq = pos.min() if not pos.empty else (end - start)

    # Internal missing time beyond one expected step; user can specify a larger freq to allow for bigger gaps, which creates negative values here that we clip to zero, i.e., no missing time counted for those.
    internal_missing = (pd.Series(ts).diff().dropna() - freq).clip(lower=pd.Timedelta(0))
    max_internal = internal_missing.max() if not internal_missing.empty else pd.Timedelta(0)

    # Leading/trailing uncovered spans (counted as-is)
    leading = max(pd.Timedelta(0), ts[0] - start)
    trailing = max(pd.Timedelta(0), end - ts[-1])

    return max(max_internal, leading, trailing)


def add_tzaware_datetime(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert 'date' + 'tz_cd' columns into a timezone-aware datetime column 'datetime_tz'.
    Example output: 2020-12-03 02:45:00+05:00
    """
    df = df.copy()
    df["date"] = pd.to_datetime(df["date"], errors="coerce")

    def localize_row(row):
        tz_name = tz_map.get(str(row["tz_cd"]).upper())
        if tz_name is None or pd.isna(row["date"]):
            return pd.NaT
        try:
            return row["date"].tz_localize(tz_name, nonexistent="shift_forward", ambiguous="NaT")
        except TypeError:
            return row["date"]  # already tz-aware
        except Exception:
            return pd.NaT

    df["date"] = df.apply(localize_row, axis=1)
    return df


def download_stations(
    param_id: str,
    station_ids: List[str],
    cache_fname: Optional[str] = None,
    datelist: pd.DatetimeIndex = pd.date_range(start='2000-01-01', end='2000-01-02'),
    stations_chunk_size: int = 10,   # download 20 sites per request
    days_per_chunk: int = 400,
) -> List[GenericObsData]:
    """
    Download IV data for many stations efficiently (station_chunk_size sites per API call).
    Uses fetch_iv_rdb() which supports multi-site requests.
    """
    station_ids = _normalize_stations(station_ids)

    def _gather_station_info(meta, site):
        """
        Lookup metadata and format it into a station_info dict.
        """
        if not meta.empty and site in meta.index:
            # Lookup metadata
            site_name = None
            lon = lat = np.nan
            if not meta.empty and site in meta.index:
                site_name = meta.loc[site].get("station_nm", None)
                lon = meta.loc[site].get("dec_long_va", np.nan)
                lat = meta.loc[site].get("dec_lat_va",  np.nan)

            # Build station info
            station_info = {
                'id': site, 'name': site_name, 'lon': lon, 'lat': lat,
                'var_name': param_id, 'var_code': param_id, 'unit': None,
            }
        else:
            raise ValueError(f"Metadata for site {site} not found.")
        return station_info

    def _standardize_df(df: pd.DataFrame) -> pd.DataFrame:
        df = add_tzaware_datetime(df)
        df = df.loc[:, ["date", "value"]].copy()
        df = df.sort_values("date").reset_index(drop=True)
        df = df.dropna(subset=["date", "value"])
        return df

    if cache_fname is not None and os.path.exists(cache_fname):
        print(f'Loading cached data from {cache_fname} ...')
        if not check_manifest(
            f"{cache_fname}.manifest",
            station_ids, param_id=param_id, start=datelist[0], end=datelist[-1]
        ):
            print("Cache manifest check failed, re-downloading data.")

        else:
            print("Cache manifest check passed, loading cached data.")

            if cache_fname.endswith((".parquet", ".pq")):
                df = load_cache_parquet(cache_fname)
                total_data = df_to_generic_obs(df)
                return total_data
            elif cache_fname.endswith(".csv"):
                df = pd.read_csv(cache_fname, parse_dates=["date"])
                # ensure tz-aware UTC
                if df["date"].dt.tz is None:
                    df["date"] = pd.to_datetime(df["date"], utc=True)
                total_data = df_to_generic_obs(df)
                return total_data

    print(f'Downloading data for {len(station_ids)} stations ...')

    start = pd.to_datetime(datelist[0]).to_pydatetime()
    end = pd.to_datetime(datelist[-1]).to_pydatetime()

    session = make_session(ignore_proxies=True)

    # Get metadata once (for names, lon/lat)
    meta = usgs_site_info(session, station_ids)
    meta = meta.set_index("site_no") if not meta.empty else pd.DataFrame()

    total_data: List[GenericObsData] = []

    # Go through one chunk of stations at a time
    missing_stations = set()
    for k, group in enumerate(chunks(station_ids, stations_chunk_size)):
        print(f"\nProcessing station group {k+1}: "
              f"{len(group) * (k) + 1} - {len(group) * (k+1)} of {len(station_ids)} ...")
        try:  # try multi-station fetch first
            df_multi = fetch_iv_rdb(session, group, param_id, start, end,
                                    days_per_chunk=days_per_chunk, timeout=(5, 180))
            if df_multi.empty:
                continue
            for site, df_site in df_multi.groupby("site_no"):
                if df_site.empty:
                    raise ValueError(f"No data for site {site} in multi-station fetch.")
                df = _standardize_df(df_site)
                station_info = _gather_station_info(meta, site)
                total_data.append(GenericObsData(station_info=station_info, df=df))
        except Exception as e:
            print(f"Error in multi-station fetch for group {k+1}: {e}")
            continue

    missing_stations = set(station_ids) - set([data.station_info['id'] for data in total_data])

    n_fallbacks_success = 0
    n_missing = len(missing_stations)
    for k, site in enumerate(missing_stations.copy()):
        print(f"\nStation {site} ({k+1} of {n_missing}) missing from multi-station fetch, "
              "attempting single-station fallback ...")
        try:
            df_single = fetch_iv_rdb(session, [site], param_id, start, end, days_per_chunk=10, timeout=(5, 180))
            if not df_single.empty:
                df = _standardize_df(df_single)
                station_info = _gather_station_info(meta, site)
                total_data.append(GenericObsData(station_info=station_info, df=df))

                n_fallbacks_success += 1
                missing_stations.remove(site)

                print(f"{site} succeeded in single-station fetch.")
            else:
                print(f"Warning: Station {site} returned no data in single-station fetch.\n")
        except Exception as e:
            print(f"Error in single-station fetch for station {site}: {e}\n")
            continue
        print(f"successful fallbacks: ({n_fallbacks_success} out of {n_missing}).\n")

    print(f'Number of stations with available data: {len(total_data)}')
    failed_stations = set(station_ids) - set([data.station_info['id'] for data in total_data])
    print(f'\n{len(failed_stations)} failed stations: {sorted(failed_stations)}')

    # Optional cache
    if cache_fname is not None:
        print(f"Saving downloaded data to {cache_fname} ...")
        ext = Path(cache_fname).suffix.lower()

        write_manifest(f"{cache_fname}.manifest", param_id, start, end, station_ids)

        if ext in (".parquet", ".pq"):
            save_cache_parquet(total_data, cache_fname)
        elif ext == ".csv":
            save_cache_csv(total_data, cache_fname)
        else:
            # fallback to parquet if no/unknown extension
            save_cache_parquet(total_data, cache_fname if ext else cache_fname + ".parquet")

    return total_data


# -------------------- Higher-level helpers (unchanged or minimal edits) --------------------

def write_time_average(input_data: List[GenericObsData],
                       param_id: Optional[str] = None,
                       unit_conv: float = 1.0,
                       outfilename: Optional[str] = None):
    """
    Writes mean of each station to 'x y z id' lines.
    If param_id == '00095' (specific conductance), convert to salinity via GSW using mS/cm at 25C, p=0 dbar.
    """
    if outfilename is None:
        raise ValueError("outfilename must be provided")
    with open(outfilename, 'w') as fout:
        for data in input_data:
            x = data.station_info.get('lon', np.nan)
            y = data.station_info.get('lat', np.nan)
            z = data.df['value'].mean()
            if param_id == '00095':  # convert conductance to salinity
                z = gsw.conversions.SP_from_C(z * unit_conv, 25, 0)
            if np.isnan(z):
                continue
            fout.write(f"{x} {y} {z} {data.station_info.get('id','')}\n")


def all_states_time_average(param_id=None, unit_conv=1.0, states=None,
                            datelist=pd.date_range(start='2000-01-01', end='2000-01-02'),
                            outfilename=None):
    if outfilename is None:
        raise ValueError("outfilename must be provided")

    if os.path.exists(outfilename):
        print(f'{outfilename} exists, skipping ...')
        return

    if states is None:
        states = ecgc_states

    # Use parameter filter to narrow site listing and speed things up
    station_info_df = get_usgs_stations_from_state(states, data_type='iv', parameter=param_id)
    station_ids = station_info_df["site_no"].tolist()

    os.makedirs(os.path.dirname(outfilename), exist_ok=True)

    total_data = download_stations(param_id=param_id,
                                   station_ids=station_ids,
                                   datelist=datelist,
                                   stations_chunk_size=50,
                                   days_per_chunk=3)
    write_time_average(input_data=total_data, param_id=param_id, unit_conv=unit_conv, outfilename=outfilename)


def convert_to_ObsData(total_data: List[GenericObsData], cache_fname: Optional[str] = None) -> ObsData:
    """
    Convert the downloaded "total_data" to ObsData/Station (for legacy scripts).
    """
    stations = []
    for data in total_data:
        my_station = Station(fname=None)
        my_station.id = data.station_info.get('id')
        my_station.lon = data.station_info.get('lon')
        my_station.lat = data.station_info.get('lat')
        my_station.df = data.df.set_index(pd.DatetimeIndex(pd.to_datetime(data.df["date"], utc=True)))
        stations.append(my_station)

    obs_data = ObsData(fname=None)
    obs_data.stations = stations
    obs_data.set_fID2id()

    if cache_fname is not None:
        obs_data.saved_file = cache_fname
        obs_data.save(obs_data.saved_file)

    return obs_data


def sample_get_stations():
    """
    Example usage of get_usgs_stations_from_state function
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
    # save station id and names in csv
    station_ids = [data.station_info['id'] for data in total_data]
    station_names = [data.station_info['name'] for data in total_data]
    # strip leading and trailing spaces and quotes
    station_ids = [sid.strip().strip('"') for sid in station_ids]
    station_names = [sname.strip().strip('"') for sname in station_names]
    # write to csv, using ';' as separator
    df = pd.DataFrame({
        'site_no': station_ids, 'station_nm': station_names
    }).drop_duplicates(subset='site_no')
    df.to_csv('station_info.csv', index=False, sep=';')
    print(df.head())


def sample_single_station():
    """
    Example usage of download_stations function for a single station
    """
    sess = make_session(ignore_proxies=True)
    df = fetch_iv_rdb(
        sess, ["01364500"], "00060",
        pd.Timestamp("2017-11-30"), pd.Timestamp("2019-01-03"),
        days_per_chunk=1, timeout=(5, 180))
    print(df.head())


def sample_read_download_params():
    """'/sciclone/schism10/feiye/Cache//usgs_stofs3d_atl_0006
    Example usage of reading usgs_download_params.json and downloading data
    """
    import json
    param_dict = json.load(open(
        '/sciclone/schism10/feiye/STOFS3D-v7.3/I18/Source_sink/USGS_adjusted_sources/Diag/usgs_download_params.json'
    ))

    data = download_stations(
        param_id=param_dict["param_id"],
        station_ids=param_dict["station_ids"],
        cache_fname='/sciclone/schism10/feiye/Cache//usgs_stofs3d_atl_00060_20191130_20210102.pq',
        datelist=pd.date_range(start="2022-11-30", end="2024-01-02")
    )
    print(f'Downloaded {len(data)} stations.')

    return data


if __name__ == "__main__":
    total_data = sample_read_download_params()
    for data in total_data:
        max_gap = detect_data_gap(
            data.df.set_index(pd.DatetimeIndex(pd.to_datetime(data.df["date"], utc=True))),
            start=pd.Timestamp("2022-12-01", tz='UTC'), end=pd.Timestamp("2024-01-01", tz='UTC')
        )
        if max_gap > pd.Timedelta('5 days'):
            print(f"Station {data.station_info['id']} max data gap: {max_gap}")
    print("Done.")
