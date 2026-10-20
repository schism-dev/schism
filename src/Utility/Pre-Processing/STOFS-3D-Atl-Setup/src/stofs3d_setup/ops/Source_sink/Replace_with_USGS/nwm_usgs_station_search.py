"""Shared NWM network preparation and USGS station search helpers."""

import os
import pickle
import time
from datetime import datetime
from glob import glob
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
from schism_py_pre_post.Utilities.util import b_in_a

from stofs3d_setup.ops.USGS.ogc.source_sink import (
    get_usgs_stations_from_state,
    usgs_var_dict,
)
from stofs3d_setup.ops.USGS.ogc.download_usgs import (
    discover_monitoring_locations,
    filter_series_by_hgrid,
)
from stofs3d_setup.utils.utils import STOFS3D_ATL_STATES


def read_nwm_data(
    data_dir=None,
    fids=None,
    date_range=pd.date_range(
        "2015-06-01 00:00:00", "2015-08-01 00:00:00", freq="3h"
    ),
):
    """Read NWM discharge values for selected feature IDs and times."""
    fids = [] if fids is None else fids
    datafiles = []
    for this_date in date_range:
        files = glob(
            f'{data_dir}/*{datetime.strftime(this_date, "%Y%m%d%H")}*'
        )
        if not files:
            raise FileNotFoundError(f"no file found for {this_date}")
        if len(files) > 1:
            raise ValueError(f"multiple files found for {this_date}")
        datafiles.append(files[0])

    stream_flow = np.full((len(datafiles), len(fids)), np.nan)
    for row, this_file in enumerate(datafiles):
        print(f"reading {this_file}\n")
        with xr.open_dataset(this_file) as dataset:
            feature_ids = dataset.feature_id.data
            positions = b_in_a(a=feature_ids, b=fids)
            if np.ma.is_masked(positions):
                raise ValueError("failed to find some segments in NWM data")
            stream_flow[row, :] = dataset.streamflow.data[positions]
    return stream_flow


class UsgsSt:
    """USGS station ID and the nearby NWM segment ID."""

    def __init__(self, st_id=None, nearby_nwm_fid=None):
        self.st_id = st_id
        self.nearby_nwm_fid = nearby_nwm_fid


class SearchPath:
    """NWM segment IDs visited while searching in one direction."""

    def __init__(self):
        self._seg_id = []
        self._order = []

    def forward_path(self, seg_id, seg_order):
        self._seg_id.append(seg_id)
        self._order.append(seg_order)

    def reset_path(self):
        self._seg_id = []
        self._order = []

    def writer(self, fname, shp_gpd=None):
        with open(fname, "w", encoding="utf-8") as output:
            output.write("id fid lon lat \n")
            for index, feature_id in enumerate(self._seg_id):
                output.write(
                    f"{index + 1} {feature_id} "
                    f"{shp_gpd.loc[feature_id].lon} "
                    f"{shp_gpd.loc[feature_id].lat}\n"
                )


class Vsource:
    """Source location and station-search state for one SCHISM source."""

    def __init__(self, xyz=None, hgrid_ie=0, nwm_fid=0, df=None):
        self.xyz = [0.0, 0.0, 0.0] if xyz is None else xyz
        self.usgs_st = []
        self.nwm_fid = nwm_fid
        self.upstream_path = SearchPath()
        self.downstream_path = SearchPath()
        self.hgrid_ie = hgrid_ie
        self.df = pd.DataFrame() if df is None else df

    @property
    def df(self):
        return self._df

    @df.setter
    def df(self, dataframe):
        if not isinstance(dataframe, pd.DataFrame):
            raise ValueError("needs a DataFrame with a 'datetime' column")
        if "datetime" in dataframe:
            self._df = dataframe.set_index(
                pd.DatetimeIndex(pd.to_datetime(dataframe["datetime"], utc=True))
            )
        else:
            self._df = dataframe


def find_usgs_along_nwm(
    iupstream=True,
    starting_seg_id=None,
    vsource=None,
    total_seg_length=0.0,
    order=0,
    nwm_shp=None,
):
    """Search up- or downstream to the first gaged segment on each branch."""
    search_distance_limit_m = 30e3
    if starting_seg_id is None or starting_seg_id not in nwm_shp.index:
        if starting_seg_id is not None:
            print(
                f"starting segment {starting_seg_id} not found in the "
                "shapefile, skipping..."
            )
        return

    path = vsource.upstream_path if iupstream else vsource.downstream_path
    path.forward_path(seg_id=starting_seg_id, seg_order=order)
    segment = nwm_shp.loc[starting_seg_id]
    total_seg_length += segment.Length
    if total_seg_length > search_distance_limit_m:
        return
    if segment.gages is not None and segment.gages != "":
        vsource.usgs_st.append(
            UsgsSt(st_id=segment.gages, nearby_nwm_fid=starting_seg_id)
        )
        return

    next_ids = segment["from"] if iupstream else [segment["to"]]
    for next_id in next_ids:
        find_usgs_along_nwm(
            iupstream=iupstream,
            starting_seg_id=next_id,
            vsource=vsource,
            total_seg_length=total_seg_length,
            order=order + 1,
            nwm_shp=nwm_shp,
        )


def prepare_usgs_stations(
    states=None,
    diag_output=None,
    parameter=None,
    start="1800-01-01T00:00:00Z",
    end="2200-01-01T00:00:00Z",
    locations=None,
):
    """Return screened USGS station IDs and coordinates."""
    if states is None:
        raise ValueError("a list of states must be provided")
    station_info = get_usgs_stations_from_state(
        states=states,
        parameter=parameter or usgs_var_dict["streamflow"]["id"],
        start=start,
        end=end,
        locations=locations,
    )
    station_ids = station_info["site_no"].to_numpy()
    station_lon = station_info["dec_long_va"].to_numpy()
    station_lat = station_info["dec_lat_va"].to_numpy()
    valid = (station_lon != "") & (station_lat != "")
    station_ids = station_ids[valid]
    station_lon = station_lon[valid].astype(float)
    station_lat = station_lat[valid].astype(float)
    _, unique_idx = np.unique(station_ids, return_index=True)
    station_ids = station_ids[unique_idx]
    station_lon = station_lon[unique_idx]
    station_lat = station_lat[unique_idx]

    if diag_output is not None:
        with open(diag_output, "w", encoding="utf-8") as output:
            output.write("lon,lat,id\n")
            for index in range(len(station_ids)):
                output.write(
                    f"{station_lon[index]},{station_lat[index]},"
                    f"{station_ids[index]}\n"
                )
    return station_ids, np.c_[station_lon, station_lat]


def prepare_nwm_usgs_station_search(
    states,
    nwm_shapefile,
    usgs_cache_folder,
    *,
    manual_feature_station_links=None,
    station_diagnostic_path=None,
    association_diagnostic_path=None,
    invalid_station_diagnostic_path="invalid_usgs_station_coordinates.csv",
    save_station_cache=True,
    parameter=None,
    start="1800-01-01T00:00:00Z",
    end="2200-01-01T00:00:00Z",
    hgrid=None,
    domain_buffer_km=100.0,
):
    """Screen domain stations by parameter, then associate them with NWM reaches."""
    states = list(states or STOFS3D_ATL_STATES)
    cache_dir = Path(usgs_cache_folder)
    cache_dir.mkdir(parents=True, exist_ok=True)
    station_cache = cache_dir / (
        f'usgs_station_inventory_{"_".join(sorted(states))}.csv'
    )
    if station_cache.is_file():
        locations = pd.read_csv(station_cache, dtype={"station_id": str})
    else:
        locations = discover_monitoring_locations(states=states)
        if save_station_cache:
            locations.to_csv(station_cache, index=False)

    locations = filter_series_by_hgrid(
        locations, hgrid, domain_buffer_km=domain_buffer_km
    )

    station_ids, station_coords = prepare_usgs_stations(
        states=states,
        parameter=parameter or usgs_var_dict["streamflow"]["id"],
        start=start,
        end=end,
        locations=locations,
        diag_output=(
            None
            if station_diagnostic_path is None
            else str(station_diagnostic_path)
        ),
    )

    nwm_shp = preprocess_nwm_shp(str(nwm_shapefile))
    nwm_shp = associate_poi_with_nwm(
        nwm_shp,
        poi=np.asarray(station_coords, dtype=float),
        poi_names=np.asarray(station_ids, dtype=str).tolist(),
        poi_label="gages",
        diag_output=(
            None
            if association_diagnostic_path is None
            else association_diagnostic_path
        ),
        invalid_poi_output=invalid_station_diagnostic_path,
    )
    for feature_id, station_id in (manual_feature_station_links or {}).items():
        matches = nwm_shp["featureID"] == int(feature_id)
        if np.any(matches):
            nwm_shp.loc[matches, "gages"] = str(station_id)
    return nwm_shp


def associate_poi_with_nwm(
    nwm_shp,
    poi: np.ndarray,
    poi_names: list | None = None,
    poi_label: str = "poi",
    diag_output: str | None = None,
    invalid_poi_output: str | None = "invalid_usgs_station_coordinates.csv",
) -> gpd.GeoDataFrame:
    """Associate each station with its nearest eligible NWM flowline."""
    poi = np.asarray(poi, dtype=float)
    if poi_names is None:
        poi_names = [str(index) for index in range(len(poi))]
    poi_names = np.asarray(poi_names)
    segment_centers = np.c_[
        nwm_shp.geometry.centroid.x,
        nwm_shp.geometry.centroid.y,
    ]
    segment_tree = cKDTree(segment_centers)
    bad_poi = ~np.isfinite(poi).all(axis=1)
    if bad_poi.any():
        print("Invalid USGS stations:")
        for name, coordinates in zip(poi_names[bad_poi], poi[bad_poi]):
            print(f"  station={name}, lon={coordinates[0]}, lat={coordinates[1]}")
        if invalid_poi_output is not None:
            pd.DataFrame(
                {
                    "station_id": poi_names[bad_poi],
                    "longitude": poi[bad_poi, 0],
                    "latitude": poi[bad_poi, 1],
                }
            ).to_csv(invalid_poi_output, index=False)
        poi = poi[~bad_poi]
        poi_names = poi_names[~bad_poi]

    print(f"poi shape: {poi.shape}")
    print(f"NaN in poi: {np.isnan(poi).sum()}")
    print(f"Inf in poi: {np.isinf(poi).sum()}")
    distances, indices = segment_tree.query(poi, k=10)
    distances[distances > 0.1] = np.inf
    indices[distances == np.inf] = -1
    poi_points = gpd.points_from_xy(poi[:, 0], poi[:, 1])
    poi_to_segment = np.full(len(poi), -1, dtype=int)
    poi_to_distance = np.full(len(poi), np.inf)

    for poi_index in range(len(poi)):
        candidates = []
        for neighbor_index in range(10):
            segment_position = indices[poi_index, neighbor_index]
            if segment_position == -1:
                continue
            segment_id = nwm_shp.index[segment_position]
            distance = nwm_shp.loc[segment_id, "geometry"].distance(
                poi_points[poi_index]
            )
            if distance <= 0.01:
                candidates.append((segment_id, distance))
        if candidates:
            segment_id, distance = min(candidates, key=lambda pair: pair[1])
            poi_to_segment[poi_index] = nwm_shp.index.get_loc(segment_id)
            poi_to_distance[poi_index] = distance

    valid = np.flatnonzero(poi_to_segment != -1)
    matches = pd.DataFrame(
        {
            "poi": [poi_names[index] for index in valid],
            "seg_pos": [poi_to_segment[index] for index in valid],
            "dist": [poi_to_distance[index] for index in valid],
        }
    )
    segment_ids = nwm_shp.index.to_numpy()
    matches["seg_idx"] = segment_ids[matches["seg_pos"].values]
    duplicates = matches[
        matches.duplicated("seg_idx", keep=False)
    ].sort_values(["seg_idx", "dist"])
    for segment_id, group in duplicates.groupby("seg_idx"):
        feature_id = nwm_shp.at[segment_id, "featureID"]
        print(f"warning: NWM segment {feature_id} matched {len(group)} POIs:")
        for _, row in group.iterrows():
            print(f"    poi={row.poi}, dist={row.dist:.6f}")
        print(
            f"    choosing the closest one: poi={group.iloc[0].poi}, "
            f"dist={group.iloc[0].dist:.6f}\n"
        )
    best_matches = (
        matches.sort_values(["seg_idx", "dist"])
        .groupby("seg_idx", as_index=False)
        .first()
    )
    nwm_shp[poi_label] = None
    nwm_shp.loc[best_matches["seg_idx"].values, poi_label] = (
        best_matches["poi"].values
    )
    print(f'{nwm_shp[poi_label].notnull().sum()} {poi_label} associated with NWM segments')

    if diag_output is not None:
        with open(diag_output, "w", encoding="utf-8") as output:
            for index, segment_position in enumerate(poi_to_segment):
                if segment_position != -1:
                    segment_id = nwm_shp.index[segment_position]
                    output.write(
                        f"{poi[index, 0]}, {poi[index, 1]}, {poi_names[index]}, "
                        f"{nwm_shp.loc[segment_id].featureID}\n"
                    )
        nwm_shp.drop(
            columns=["featureID", "from", "to"], errors="ignore"
        ).to_file(f"{diag_output}.gpkg", layer="nwm_with_poi", driver="GPKG")
    return nwm_shp


def preprocess_nwm_shp(f_shapefile):
    """Load the NWM network and attach upstream segment IDs."""
    started_at = time.time()
    cache_path = f"{f_shapefile}.pkl"
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as cache_file:
            nwm_shp = pickle.load(cache_file)
    else:
        nwm_shp = gpd.read_file(f_shapefile)
        index_by_feature_id = {
            record.featureID: index
            for index, record in nwm_shp.iterrows()
        }
        upstream_segments = [[] for _ in range(len(nwm_shp))]
        for _, record in nwm_shp.iterrows():
            try:
                upstream_segments[
                    index_by_feature_id[record.to]
                ].append(record.featureID)
            except KeyError:
                continue
        nwm_shp["from"] = upstream_segments
        nwm_shp = nwm_shp.set_index("featureID")
        nwm_shp["featureID"] = nwm_shp.index
        with open(cache_path, "wb") as cache_file:
            pickle.dump(nwm_shp, cache_file)
    print(
        "---------------loading and processing shapefile took: "
        f"{time.time() - started_at} s ---------------\n"
    )
    return nwm_shp
