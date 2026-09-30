"""Automatic USGS temperature corrections for SCHISM sources."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from stofs3d_setup.ops.Source_sink.Replace_with_USGS.download_usgs import (
    download_stations,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.replace_with_obs import (
    Vsource,
    associate_poi_with_nwm,
    find_usgs_along_nwm,
    prepare_usgs_stations,
    preprocess_nwm_shp,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    MANUAL_NWM_TO_USGS_FLOW,
    USGS_TEMPERATURE_PARAMETER_ID,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.usgs_series import (
    _as_utc_timestamp,
    _extract_usgs_values,
    _interpolate_usgs_with_original_fallback,
    _model_datetimes,
    _station_id_from_record,
)
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _check_matching_time,
)
from stofs3d_setup.utils.utils import STOFS3D_ATL_STATES


def _load_relocated_source_fids(
    relocated_source_sink_dir: str | Path,
) -> dict[int, list[int]]:
    """Read relocated element-to-NWM-feature mapping from sources.json."""
    mapping_file = Path(relocated_source_sink_dir) / "sources.json"
    if not mapping_file.is_file():
        raise FileNotFoundError(
            "Relocated sources.json is required for automatic temperature "
            f"replacement: {mapping_file}"
        )

    with mapping_file.open("r", encoding="utf-8") as f:
        raw = json.load(f)

    mapping: dict[int, list[int]] = {}
    for ele, fids in raw.items():
        if isinstance(fids, (int, np.integer, str)):
            fids = [fids]
        mapping[int(ele)] = [int(fid) for fid in fids]

    return mapping

def _prepare_nwm_usgs_station_search(
    states,
    nwm_shapefile,
    usgs_cache_folder,
    diagnostics_dir: Path | None,
):
    """
    Build the NWM network with USGS streamflow stations attached.

    This uses the same station inventory, NWM association, and upstream search
    functions as source_nwm2usgs().
    """
    states = list(states or STOFS3D_ATL_STATES)
    cache_dir = Path(usgs_cache_folder)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cache_name = f'usgs_states_{"_".join(np.sort(states))}.txt'
    station_cache = cache_dir / cache_name

    if station_cache.is_file():
        station_df = pd.read_csv(station_cache, dtype=object)
        station_ids = station_df["id"].astype(str).to_numpy()
        station_coords = (
            station_df[["lon", "lat"]].to_numpy(dtype=float)
        )
    else:
        diag_output = (
            None
            if diagnostics_dir is None
            else str(diagnostics_dir / "all_source_usgs_stations.csv")
        )
        station_ids, station_coords = prepare_usgs_stations(
            states=states,
            diag_output=diag_output,
        )
        station_df = pd.DataFrame(
            {
                "id": np.asarray(station_ids, dtype=str),
                "lon": station_coords[:, 0],
                "lat": station_coords[:, 1],
            }
        )
        station_df.to_csv(station_cache, index=False)

    nwm_shp = preprocess_nwm_shp(str(nwm_shapefile))

    association_diag = (
        None
        if diagnostics_dir is None
        else str(diagnostics_dir / "all_source_auto_nearby_gages.txt")
    )
    nwm_shp = associate_poi_with_nwm(
        nwm_shp,
        poi=np.asarray(station_coords, dtype=float),
        poi_names=np.asarray(station_ids, dtype=str).tolist(),
        poi_label="gages",
        diag_output=association_diag,
    )

    # Retain the manual feature-to-station associations used by the flow
    # replacement workflow, including current artificial-island overrides.
    for feature_id, station_id in MANUAL_NWM_TO_USGS_FLOW.items():
        idx = nwm_shp["featureID"] == int(feature_id)
        if np.any(idx):
            nwm_shp.loc[idx, "gages"] = str(station_id)

    return nwm_shp


def _find_station_ids_for_relocated_sources(
    source_eles: list[int],
    relocated_ele_to_fids: dict[int, list[int]],
    source_time_and_data,
    xctr: np.ndarray,
    yctr: np.ndarray,
    nwm_shp,
) -> dict[int, list[str]]:
    """
    Find upstream USGS stations for every relocated source.

    Multiple NWM feature IDs assigned to one relocated element are searched.
    Station IDs are de-duplicated while preserving search order.
    """
    if source_time_and_data is None:
        return {}

    _, source_data = source_time_and_data
    source_station_ids: dict[int, list[str]] = {}

    for source_idx, source_ele in enumerate(source_eles):
        fids = relocated_ele_to_fids.get(int(source_ele), [])
        if not fids:
            source_station_ids[int(source_ele)] = []
            continue

        found: list[str] = []
        for fid in fids:
            vsource = Vsource(
                xyz=[
                    float(xctr[source_ele - 1]),
                    float(yctr[source_ele - 1]),
                    0.0,
                ],
                hgrid_ie=int(source_ele),
                nwm_fid=int(fid),
                df=pd.DataFrame(
                    {
                        "datetime": pd.date_range(
                            "2000-01-01",
                            periods=source_data.shape[0],
                            freq="1h",
                            tz="UTC",
                        ),
                        "Data": source_data[:, source_idx],
                    }
                ),
            )

            find_usgs_along_nwm(
                iupstream=True,
                starting_seg_id=int(fid),
                vsource=vsource,
                total_seg_length=0.0,
                order=0,
                nwm_shp=nwm_shp,
            )

            for station in vsource.usgs_st:
                station_id = str(station.st_id)
                if station_id not in found:
                    found.append(station_id)

        source_station_ids[int(source_ele)] = found

    return source_station_ids



def _download_temperature_station_map(
    station_ids: list[str],
    start_time,
    end_time,
    usgs_cache_folder,
) -> dict[str, pd.Series]:
    """Original bulk temperature download; no custom chunking."""
    station_ids = sorted({str(station_id) for station_id in station_ids})
    if not station_ids:
        return {}

    cache_dir = Path(usgs_cache_folder)
    cache_dir.mkdir(parents=True, exist_ok=True)

    start = _as_utc_timestamp(start_time)
    end = _as_utc_timestamp(end_time)
    padded_start = start - pd.Timedelta(days=1)
    padded_end = end + pd.Timedelta(days=1)

    cache_file = cache_dir / (
        "all_source_usgs_temperature_00010_"
        f"{padded_start.strftime('%Y%m%d')}_"
        f"{padded_end.strftime('%Y%m%d')}.pq"
    )

    try:
        records = download_stations(
            param_id=USGS_TEMPERATURE_PARAMETER_ID,
            station_ids=station_ids,
            cache_fname=str(cache_file),
            datelist=pd.date_range(
                start=padded_start.tz_localize(None),
                end=padded_end.tz_localize(None),
            ),
        )
    except Exception as exc:
        print(
            "[ARTIFICIAL ISLAND PATCH] warning: bulk USGS temperature "
            f"download failed: {exc}"
        )
        return {}

    result = {}
    for record in records:
        station_id = _station_id_from_record(record)
        if not station_id:
            continue
        try:
            series = _extract_usgs_values(record)
        except Exception as exc:
            print(
                "[ARTIFICIAL ISLAND PATCH] warning: could not parse "
                f"temperature for station {station_id}: {exc}"
            )
            continue
        if not series.empty:
            result[station_id] = series

    return result


def _replace_all_relocated_source_temperatures(
    source_eles: list[int],
    source_time_and_data,
    msource_data_list,
    xctr: np.ndarray,
    yctr: np.ndarray,
    relocated_source_sink_dir,
    start_time,
    usgs_cache_folder,
    nwm_shapefile,
    states,
    diagnostics_dir: Path | None = None,
):
    """
    Automatically replace temperature for all relocated sources.

    For each relocated source:
      relocated element -> NWM feature IDs from relocated sources.json
      -> upstream USGS stations using source_nwm2usgs search logic
      -> first station with usable 00010 data
      -> blend with existing temperature; long gaps retain existing values.

    Flow is not changed in this step.
    """
    if source_time_and_data is None or not source_eles:
        return msource_data_list, 0

    if not msource_data_list:
        print(
            "[ARTIFICIAL ISLAND PATCH] warning: no msource tracer exists; "
            "automatic temperature replacement skipped."
        )
        return msource_data_list, 0

    source_time, source_data = source_time_and_data
    tracer_time, temperature_data = msource_data_list[0]
    _check_matching_time(
        source_time,
        tracer_time,
        "automatic all-source temperature replacement",
    )

    relocated_ele_to_fids = _load_relocated_source_fids(
        relocated_source_sink_dir
    )
    nwm_shp = _prepare_nwm_usgs_station_search(
        states=states,
        nwm_shapefile=nwm_shapefile,
        usgs_cache_folder=usgs_cache_folder,
        diagnostics_dir=diagnostics_dir,
    )
    source_station_ids = _find_station_ids_for_relocated_sources(
        source_eles=source_eles,
        relocated_ele_to_fids=relocated_ele_to_fids,
        source_time_and_data=source_time_and_data,
        xctr=xctr,
        yctr=yctr,
        nwm_shp=nwm_shp,
    )

    all_station_ids = [
        station_id
        for ids in source_station_ids.values()
        for station_id in ids
    ]

    target_datetime = _model_datetimes(start_time, source_time)
    temperature_station_map = _download_temperature_station_map(
        station_ids=all_station_ids,
        start_time=target_datetime[0],
        end_time=target_datetime[-1],
        usgs_cache_folder=usgs_cache_folder,
    )

    replaced_count = 0
    total_usgs_records = 0

    for source_idx, source_ele in enumerate(source_eles):
        station_ids = source_station_ids.get(int(source_ele), [])
        used_station = None
        used_mask = None

        for station_id in station_ids:
            series = temperature_station_map.get(str(station_id))
            if series is None:
                continue

            existing_temperature = temperature_data[:, source_idx].copy()
            replaced_temperature, use_usgs = (
                _interpolate_usgs_with_original_fallback(
                    series=series,
                    target_time=target_datetime,
                    original_values=existing_temperature,
                    scale=1.0,
                )
            )

            if not np.any(use_usgs):
                continue

            temperature_data[:, source_idx] = replaced_temperature
            used_station = str(station_id)
            used_mask = use_usgs
            break

        if used_station is not None:
            replaced_count += 1
            total_usgs_records += int(used_mask.sum())
            print(
                "[ARTIFICIAL ISLAND PATCH] automatic temperature: "
                f"source {source_idx}_{source_ele} using USGS "
                f"{used_station} for {int(used_mask.sum())}/"
                f"{len(used_mask)} records."
            )

    msource_data_list[0] = (tracer_time, temperature_data)

    print(
        "[ARTIFICIAL ISLAND PATCH] automatic temperature replacement "
        f"completed for {replaced_count}/{len(source_eles)} relocated "
        f"source(s), {total_usgs_records} USGS-backed records total."
    )

    return msource_data_list, replaced_count
