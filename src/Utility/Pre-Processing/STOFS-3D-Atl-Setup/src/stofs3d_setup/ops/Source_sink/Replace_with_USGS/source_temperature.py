"""Automatic USGS temperature corrections for SCHISM sources."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd

from stofs3d_setup.ops.USGS.ogc.source_sink import (
    RequiredStationDownloadError,
    download_stations,
)
from stofs3d_setup.ops.USGS.ogc.download_policy import (
    USGS_MAX_MISSING_GAP_DAYS,
)
from stofs3d_setup.ops.USGS.ogc.download_usgs import MIN_WATER_TEMPERATURE_C
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.nwm_usgs_station_search import (
    Vsource,
    find_usgs_along_nwm,
    prepare_nwm_usgs_station_search,
    read_nwm_data,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    NWM_TO_USGS_TEMPERATURE_SEARCH,
    OPTIONAL_TEMPERATURE_STATIONS,
    REQUIRED_TEMPERATURE_STATIONS,
    TEMPERATURE_SEARCH_DIRECTIONS,
    USGS_TEMPERATURE_PARAMETER_ID,
)
from stofs3d_setup.ops.USGS.ogc.series import (
    _as_utc_timestamp,
    _extract_usgs_values,
    _interpolate_usgs_with_original_fallback,
    _model_datetimes,
    _station_id_from_record,
)
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _build_source_sink,
    _check_matching_time,
    _copy_sink_components,
    _copy_source_components,
)
from stofs3d_setup.ops.Source_sink.spatial_selection import (
    _compute_grid_centers,
)
from stofs3d_setup.utils.utils import STOFS3D_ATL_STATES


FIRST_USABLE = "first_usable"
DISCHARGE_WEIGHTED = "discharge_weighted"
TEMPERATURE_POOLING_MODES = {FIRST_USABLE, DISCHARGE_WEIGHTED}
AMBIENT_TEMPERATURE_SENTINEL = -9999.0
MAX_TEMPERATURE_MISSING_GAP = pd.Timedelta(
    days=USGS_MAX_MISSING_GAP_DAYS
)


@dataclass(frozen=True)
class TemperatureStationCandidate:
    """USGS temperature station and its associated NWM flow segment."""

    station_id: str
    nwm_feature_id: int


@dataclass(frozen=True)
class TemperatureReplacementInputs:
    """Shared observations and topology used by the pooling methods."""

    target_time: pd.DatetimeIndex
    candidates_by_source: dict[int, list[TemperatureStationCandidate]]
    temperature_by_station: dict[str, pd.Series]
    discharge_by_feature: dict[int, np.ndarray]


def _load_source_feature_mapping(
    source_mapping_dir: str | Path,
) -> dict[int, list[int]]:
    """Read source-element-to-NWM-feature mapping from sources.json."""
    mapping_file = Path(source_mapping_dir) / "sources.json"
    if not mapping_file.is_file():
        raise FileNotFoundError(
            "sources.json is required for automatic temperature "
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
    start_time,
    end_time,
    hgrid,
):
    """Build the shared NWM/USGS station network for temperature search."""
    diagnostics_dir = (
        None if diagnostics_dir is None else Path(diagnostics_dir)
    )
    return prepare_nwm_usgs_station_search(
        states=states,
        nwm_shapefile=nwm_shapefile,
        usgs_cache_folder=usgs_cache_folder,
        parameter=USGS_TEMPERATURE_PARAMETER_ID,
        start=start_time,
        end=end_time,
        hgrid=hgrid,
        manual_feature_station_links=NWM_TO_USGS_TEMPERATURE_SEARCH,
        station_diagnostic_path=(
            None
            if diagnostics_dir is None
            else diagnostics_dir / "all_source_usgs_stations.csv"
        ),
        association_diagnostic_path=(
            None
            if diagnostics_dir is None
            else diagnostics_dir / "all_source_auto_nearby_gages.txt"
        ),
        invalid_station_diagnostic_path=(
            None
            if diagnostics_dir is None
            else diagnostics_dir / "invalid_usgs_station_coordinates.csv"
        ),
        save_station_cache=True,
    )


def _find_temperature_station_candidates(
    source_eles: list[int],
    source_element_to_fids: dict[int, list[int]],
    source_time_and_data,
    xctr: np.ndarray,
    yctr: np.ndarray,
    nwm_shp,
    search_directions: tuple[str, ...] | list[str] = (
        TEMPERATURE_SEARCH_DIRECTIONS
    ),
) -> dict[int, list[TemperatureStationCandidate]]:
    """
    Find USGS stations in configured NWM network directions for each source.

    Multiple NWM feature IDs assigned to one source element are searched.
    Manually linked stations for those features are tried first; other station
    IDs are de-duplicated while preserving configured direction and network
    search order.
    """
    invalid_directions = set(search_directions) - {"upstream", "downstream"}
    if not search_directions or invalid_directions:
        raise ValueError(
            "Temperature search directions must include upstream and/or "
            f"downstream; received {search_directions!r}"
        )

    if source_time_and_data is None:
        return {}

    _, source_data = source_time_and_data
    source_candidates: dict[int, list[TemperatureStationCandidate]] = {}

    for source_idx, source_ele in enumerate(source_eles):
        fids = source_element_to_fids.get(int(source_ele), [])
        if not fids:
            source_candidates[int(source_ele)] = []
            continue

        preferred_station_ids = {
            str(NWM_TO_USGS_TEMPERATURE_SEARCH[int(fid)])
            for fid in fids
            if int(fid) in NWM_TO_USGS_TEMPERATURE_SEARCH
        }

        found: list[TemperatureStationCandidate] = []
        found_station_ids: set[str] = set()
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

            for direction in search_directions:
                find_usgs_along_nwm(
                    iupstream=direction == "upstream",
                    starting_seg_id=int(fid),
                    vsource=vsource,
                    total_seg_length=0.0,
                    order=0,
                    nwm_shp=nwm_shp,
                )

            for station in vsource.usgs_st:
                station_id = str(station.st_id)
                if station_id in found_station_ids:
                    continue
                found.append(
                    TemperatureStationCandidate(
                        station_id=station_id,
                        nwm_feature_id=int(station.nearby_nwm_fid),
                    )
                )
                found_station_ids.add(station_id)

        found.sort(
            key=lambda candidate: (
                candidate.station_id not in preferred_station_ids
            )
        )
        source_candidates[int(source_ele)] = found

    return source_candidates


def _download_temperature_station_map(
    station_ids: list[str],
    start_time,
    end_time,
    usgs_cache_folder,
    required_station_ids: list[str] | None = None,
    optional_station_ids: list[str] | None = None,
    warning_log_path: str | Path | None = None,
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
            required_station_ids=required_station_ids,
            optional_station_ids=optional_station_ids,
            validation_start=start,
            validation_end=end,
            warning_log_path=warning_log_path,
            warning_input_files=(
                "Source_sink/msource.th",
                "Source_sink/source.nc",
            ),
            warning_forcing_name="temperature",
        )
    except RequiredStationDownloadError:
        raise
    except Exception as exc:
        print(
            "[SOURCE TEMPERATURE] warning: bulk USGS temperature "
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
                "[SOURCE TEMPERATURE] warning: could not parse "
                f"temperature for station {station_id}: {exc}"
            )
            continue
        if not series.empty:
            result[station_id] = series

    return result


def _validate_temperature_pooling(pooling: str) -> str:
    """Return a normalized supported temperature-pooling mode."""
    pooling = str(pooling).strip().lower()
    if pooling not in TEMPERATURE_POOLING_MODES:
        raise ValueError(
            "source_temperature_pooling must be one of "
            f"{sorted(TEMPERATURE_POOLING_MODES)}, got {pooling!r}"
        )
    return pooling


def _is_ambient_temperature(values: np.ndarray) -> np.ndarray:
    """Identify the SCHISM ambient-temperature sentinel."""
    return np.isclose(
        np.asarray(values, dtype=float),
        AMBIENT_TEMPERATURE_SENTINEL,
        rtol=0.0,
        atol=1.0e-8,
    )


def _mixed_temperature_columns(temperature_data: np.ndarray) -> np.ndarray:
    """Return columns that are neither all ambient nor all finite numeric."""
    temperature_data = np.asarray(temperature_data, dtype=float)
    if temperature_data.ndim == 1:
        temperature_data = temperature_data.reshape(-1, 1)
    ambient = _is_ambient_temperature(temperature_data)
    all_ambient = np.all(ambient, axis=0)
    all_numeric = np.all(~ambient & np.isfinite(temperature_data), axis=0)
    return ~(all_ambient | all_numeric)


def _blend_complete_temperature_or_preserve(
    series: pd.Series | None,
    target_time: pd.DatetimeIndex,
    original_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Apply a station series only when the final column is temporally safe."""
    original_values = np.asarray(original_values, dtype=float).reshape(-1)
    if _mixed_temperature_columns(original_values)[0]:
        raise ValueError(
            "discharge_weighted temperature mode received a source column "
            "that is not uniformly all-ambient or all-finite-numeric"
        )

    blended = (
        None
        if series is None
        else _complete_temperature_series(series, target_time)
    )
    if blended is None:
        return (
            original_values.copy(),
            np.zeros(len(original_values), dtype=bool),
            False,
        )
    return blended, np.ones(len(original_values), dtype=bool), True


def _read_nwm_discharge_map(
    candidates_by_source: dict[int, list[TemperatureStationCandidate]],
    target_time: pd.DatetimeIndex,
    nwm_data_dir: str | Path,
) -> dict[int, np.ndarray]:
    """Read and interpolate NWM flow for all station-associated segments."""
    feature_ids = sorted(
        {
            candidate.nwm_feature_id
            for candidates in candidates_by_source.values()
            for candidate in candidates
        }
    )
    if not feature_ids:
        return {}
    if nwm_data_dir is None:
        raise ValueError(
            "nwm_data_dir is required for discharge_weighted temperature"
        )

    target_time = pd.DatetimeIndex(target_time)
    target_time_naive = (
        target_time.tz_convert(None)
        if target_time.tz is not None
        else target_time
    )
    nwm_time = pd.date_range(
        start=target_time_naive[0],
        end=target_time_naive[-1],
        freq="3h",
    )
    nwm_flow = read_nwm_data(
        data_dir=str(nwm_data_dir),
        fids=feature_ids,
        date_range=nwm_time,
    )

    target_seconds = (
        target_time_naive - target_time_naive[0]
    ).total_seconds().to_numpy(dtype=float)
    nwm_seconds = (
        nwm_time - target_time_naive[0]
    ).total_seconds().to_numpy(dtype=float)

    discharge_by_feature: dict[int, np.ndarray] = {}
    for feature_idx, feature_id in enumerate(feature_ids):
        flow = np.asarray(nwm_flow[:, feature_idx], dtype=float)
        valid = np.isfinite(flow)
        if np.count_nonzero(valid) < 2:
            discharge_by_feature[feature_id] = np.full(
                len(target_time), np.nan, dtype=float
            )
            continue
        discharge_by_feature[feature_id] = np.interp(
            target_seconds,
            nwm_seconds[valid],
            flow[valid],
            left=np.nan,
            right=np.nan,
        )
    return discharge_by_feature


def _temperature_on_model_grid(
    series: pd.Series,
    target_time: pd.DatetimeIndex,
) -> np.ndarray:
    """Return a complete station series or unsupported model times."""
    values = _complete_temperature_series(series, target_time)
    return (
        values
        if values is not None
        else np.full(len(target_time), np.nan, dtype=float)
    )


def _fill_short_composite_gaps(
    values: np.ndarray,
    target_time: pd.DatetimeIndex,
) -> np.ndarray:
    """Fill at most seven missing days in a weighted composite series."""
    values = np.asarray(values, dtype=float).reshape(-1)
    valid = np.isfinite(values)
    if not np.any(valid):
        return values
    series = pd.Series(values[valid], index=target_time[valid])
    filled = _complete_temperature_series(series, target_time)
    return values if filled is None else filled


def _discharge_weighted_temperature(
    candidates: list[TemperatureStationCandidate],
    temperature_by_station: dict[str, pd.Series],
    discharge_by_feature: dict[int, np.ndarray],
    target_time: pd.DatetimeIndex,
) -> tuple[np.ndarray | None, dict[str, float | int]]:
    """Return a complete discharge-weighted series or ``None``."""
    n_times = len(target_time)
    weighted_temperature = np.zeros(n_times, dtype=float)
    observed_discharge = np.zeros(n_times, dtype=float)
    total_discharge = np.zeros(n_times, dtype=float)
    available_station_count = np.zeros(n_times, dtype=int)

    for candidate in candidates:
        discharge = discharge_by_feature.get(candidate.nwm_feature_id)
        if discharge is None:
            continue
        discharge = np.asarray(discharge, dtype=float).reshape(-1)
        if len(discharge) != n_times:
            raise ValueError(
                "NWM discharge length does not match the model time array"
            )

        valid_discharge = np.isfinite(discharge) & (discharge > 0.0)
        total_discharge[valid_discharge] += discharge[valid_discharge]

        temperature_series = temperature_by_station.get(
            candidate.station_id
        )
        if temperature_series is None:
            continue
        temperature = _temperature_on_model_grid(
            temperature_series,
            target_time,
        )
        valid_pair = valid_discharge & np.isfinite(temperature)
        weighted_temperature[valid_pair] += (
            discharge[valid_pair] * temperature[valid_pair]
        )
        observed_discharge[valid_pair] += discharge[valid_pair]
        available_station_count[valid_pair] += 1

    composite = np.full(n_times, np.nan, dtype=float)
    usable = observed_discharge > 0.0
    composite[usable] = (
        weighted_temperature[usable] / observed_discharge[usable]
    )
    composite = _fill_short_composite_gaps(composite, target_time)

    coverage = np.full(n_times, np.nan, dtype=float)
    positive_total = total_discharge > 0.0
    coverage[positive_total] = (
        observed_discharge[positive_total]
        / total_discharge[positive_total]
    )
    finite_coverage = coverage[np.isfinite(coverage)]
    diagnostics: dict[str, float | int] = {
        "candidate_station_count": len(candidates),
        "supported_records_before_composite_fill": int(np.count_nonzero(usable)),
        "single_station_records": int(
            np.count_nonzero(available_station_count == 1)
        ),
        "minimum_discharge_coverage": (
            float(np.min(finite_coverage))
            if finite_coverage.size
            else np.nan
        ),
        "mean_discharge_coverage": (
            float(np.mean(finite_coverage))
            if finite_coverage.size
            else np.nan
        ),
        "median_discharge_coverage": (
            float(np.median(finite_coverage))
            if finite_coverage.size
            else np.nan
        ),
    }

    if not np.all(np.isfinite(composite)):
        diagnostics["complete"] = 0
        diagnostics["unsupported_records_after_fill"] = int(
            np.count_nonzero(~np.isfinite(composite))
        )
        return None, diagnostics

    diagnostics["complete"] = 1
    diagnostics["unsupported_records_after_fill"] = 0
    return composite, diagnostics


def _write_temperature_candidate_diagnostics(
    candidates_by_source: dict[int, list[TemperatureStationCandidate]],
    diagnostics_dir: Path,
) -> None:
    """Write the ordered station candidates considered for each source."""
    rows = [
        {
            "source_element": source_ele,
            "candidate_order": candidate_order,
            "station_id": candidate.station_id,
            "nwm_feature_id": candidate.nwm_feature_id,
        }
        for source_ele, candidates in candidates_by_source.items()
        for candidate_order, candidate in enumerate(candidates, start=1)
    ]
    pd.DataFrame(
        rows,
        columns=[
            "source_element",
            "candidate_order",
            "station_id",
            "nwm_feature_id",
        ],
    ).to_csv(
        diagnostics_dir / "source_temperature_candidates.csv",
        index=False,
    )


def _prepare_temperature_replacement_inputs(
    source_eles: list[int],
    source_time_and_data,
    xctr: np.ndarray,
    yctr: np.ndarray,
    source_mapping_dir,
    start_time,
    usgs_cache_folder,
    nwm_shapefile,
    states,
    diagnostics_dir: Path | None = None,
    pooling: str = FIRST_USABLE,
    nwm_data_dir: str | Path | None = None,
    warning_log_path: str | Path | None = None,
    search_directions: tuple[str, ...] | list[str] = (
        TEMPERATURE_SEARCH_DIRECTIONS
    ),
    hgrid=None,
) -> TemperatureReplacementInputs:
    """Discover stations and load time series shared by pooling methods."""
    source_time, _ = source_time_and_data
    target_time = _model_datetimes(start_time, source_time)
    source_element_to_fids = _load_source_feature_mapping(source_mapping_dir)
    nwm_shp = _prepare_nwm_usgs_station_search(
        states=states,
        nwm_shapefile=nwm_shapefile,
        usgs_cache_folder=usgs_cache_folder,
        diagnostics_dir=diagnostics_dir,
        start_time=target_time[0],
        end_time=target_time[-1],
        hgrid=hgrid,
    )
    candidates_by_source = _find_temperature_station_candidates(
        source_eles=source_eles,
        source_element_to_fids=source_element_to_fids,
        source_time_and_data=source_time_and_data,
        xctr=xctr,
        yctr=yctr,
        nwm_shp=nwm_shp,
        search_directions=search_directions,
    )
    all_station_ids = [
        candidate.station_id
        for candidates in candidates_by_source.values()
        for candidate in candidates
    ]
    temperature_by_station = _download_temperature_station_map(
        station_ids=all_station_ids,
        start_time=target_time[0],
        end_time=target_time[-1],
        usgs_cache_folder=usgs_cache_folder,
        required_station_ids=sorted(
            REQUIRED_TEMPERATURE_STATIONS
            & set(all_station_ids)
        ),
        optional_station_ids=sorted(
            OPTIONAL_TEMPERATURE_STATIONS & set(all_station_ids)
        ),
        warning_log_path=warning_log_path,
    )
    discharge_by_feature = (
        _read_nwm_discharge_map(
            candidates_by_source=candidates_by_source,
            target_time=target_time,
            nwm_data_dir=nwm_data_dir,
        )
        if pooling == DISCHARGE_WEIGHTED
        else {}
    )
    if diagnostics_dir is not None:
        _write_temperature_candidate_diagnostics(
            candidates_by_source,
            diagnostics_dir,
        )

    return TemperatureReplacementInputs(
        target_time=target_time,
        candidates_by_source=candidates_by_source,
        temperature_by_station=temperature_by_station,
        discharge_by_feature=discharge_by_feature,
    )


def _first_usable_temperature(
    candidates: list[TemperatureStationCandidate],
    temperature_by_station: dict[str, pd.Series],
    target_time: pd.DatetimeIndex,
) -> tuple[np.ndarray | None, str | None, int]:
    """Return the first candidate covering every model time."""
    for candidate in candidates:
        series = temperature_by_station.get(candidate.station_id)
        if series is None:
            continue

        replacement = _complete_temperature_series(series, target_time)
        if replacement is not None:
            return replacement, candidate.station_id, len(replacement)

    return None, None, 0


def _complete_temperature_series(
    series: pd.Series,
    target_time: pd.DatetimeIndex,
) -> np.ndarray | None:
    """Fill gaps of at most seven missing days, or reject the whole series."""
    observations = series.replace([np.inf, -np.inf], np.nan).dropna()
    observations = observations[observations >= MIN_WATER_TEMPERATURE_C]
    observations = observations[~observations.index.duplicated(keep="first")]
    observations = observations.sort_index()
    if len(observations) < 2:
        return None

    if observations.index.tz is None:
        observations.index = observations.index.tz_localize("UTC")
    else:
        observations.index = observations.index.tz_convert("UTC")
    model_time = pd.DatetimeIndex(target_time)
    if model_time.tz is None:
        model_time = model_time.tz_localize("UTC")
    else:
        model_time = model_time.tz_convert("UTC")

    values, _ = _interpolate_usgs_with_original_fallback(
        series=observations,
        target_time=model_time,
        original_values=np.full(len(model_time), np.nan),
        max_interp_gap=(
            MAX_TEMPERATURE_MISSING_GAP + pd.Timedelta(hours=1)
        ),
    )

    leading = model_time < observations.index[0]
    if leading.any():
        if observations.index[0] - model_time[0] > MAX_TEMPERATURE_MISSING_GAP:
            return None
        values[leading] = float(observations.iloc[0])

    trailing = model_time > observations.index[-1]
    if trailing.any():
        if model_time[-1] - observations.index[-1] > MAX_TEMPERATURE_MISSING_GAP:
            return None
        values[trailing] = float(observations.iloc[-1])

    return values if np.all(np.isfinite(values)) else None


def _apply_temperature_replacements(
    source_eles: list[int],
    temperature_data: np.ndarray,
    inputs: TemperatureReplacementInputs,
    pooling: str,
    diagnostics_dir: Path | None,
) -> tuple[int, int]:
    """Apply a pooling method to eligible temperature columns in place."""

    replaced_count = 0
    total_usgs_records = 0
    weighted_diagnostic_rows = []

    for source_idx, source_ele in enumerate(source_eles):
        candidates = inputs.candidates_by_source.get(int(source_ele), [])

        if pooling == FIRST_USABLE:
            replacement, used_station, used_records = (
                _first_usable_temperature(
                    candidates=candidates,
                    temperature_by_station=inputs.temperature_by_station,
                    target_time=inputs.target_time,
                )
            )
            if used_station is not None:
                temperature_data[:, source_idx] = replacement
                replaced_count += 1
                total_usgs_records += used_records
                print(
                    "[SOURCE TEMPERATURE] automatic first_usable: "
                    f"source {source_idx}_{source_ele} using USGS "
                    f"{used_station} for {used_records}/"
                    f"{len(inputs.target_time)} records."
                )
            continue

        weighted_temperature, weighted_diagnostics = (
            _discharge_weighted_temperature(
                candidates=candidates,
                temperature_by_station=inputs.temperature_by_station,
                discharge_by_feature=inputs.discharge_by_feature,
                target_time=inputs.target_time,
            )
        )
        weighted_diagnostic_rows.append(
            {
                "source_index": source_idx,
                "source_element": source_ele,
                "station_ids": ";".join(
                    candidate.station_id for candidate in candidates
                ),
                "nwm_feature_ids": ";".join(
                    str(candidate.nwm_feature_id)
                    for candidate in candidates
                ),
                **weighted_diagnostics,
            }
        )
        if weighted_temperature is None:
            print(
                "[SOURCE TEMPERATURE] discharge_weighted: retained "
                f"ambient temperature for source {source_idx}_{source_ele}; "
                "a complete numeric series could not be formed."
            )
            continue

        temperature_data[:, source_idx] = weighted_temperature
        replaced_count += 1
        total_usgs_records += len(weighted_temperature)
        print(
            "[SOURCE TEMPERATURE] discharge_weighted: replaced source "
            f"{source_idx}_{source_ele} for all "
            f"{len(weighted_temperature)} records using "
            f"{len(candidates)} candidate station(s)."
        )

    mixed_columns = np.flatnonzero(_mixed_temperature_columns(temperature_data))
    if len(mixed_columns) > 0:
        mixed_elements = [source_eles[idx] for idx in mixed_columns]
        raise ValueError(
            "Temperature replacement produced columns that are not uniformly "
            "all-ambient or all-finite-numeric for source elements "
            f"{mixed_elements}"
        )

    if pooling == DISCHARGE_WEIGHTED:
        if diagnostics_dir is not None:
            pd.DataFrame(weighted_diagnostic_rows).to_csv(
                diagnostics_dir
                / "discharge_weighted_temperature_summary.csv",
                index=False,
            )

    return replaced_count, total_usgs_records


def _apply_automatic_temperature_replacements(
    source_eles: list[int],
    source_time_and_data,
    msource_data_list,
    xctr: np.ndarray,
    yctr: np.ndarray,
    source_mapping_dir,
    start_time,
    usgs_cache_folder,
    nwm_shapefile,
    states,
    diagnostics_dir: Path | None = None,
    pooling: str = FIRST_USABLE,
    nwm_data_dir: str | Path | None = None,
    warning_log_path: str | Path | None = None,
    search_directions: tuple[str, ...] | list[str] = (
        TEMPERATURE_SEARCH_DIRECTIONS
    ),
    hgrid=None,
):
    """Apply automatic USGS temperature corrections to eligible sources.

    ``first_usable`` uses the first station with complete temperature coverage.
    ``discharge_weighted`` changes a source only when a complete numeric
    temperature series can be formed. Flow is not changed by either mode.
    """
    pooling = _validate_temperature_pooling(pooling)
    if source_time_and_data is None or not source_eles:
        return msource_data_list, 0

    if not msource_data_list:
        print(
            "[SOURCE TEMPERATURE] warning: no msource tracer exists; "
            "automatic temperature replacement skipped."
        )
        return msource_data_list, 0

    source_time, _ = source_time_and_data
    tracer_time, temperature_data = msource_data_list[0]
    _check_matching_time(
        source_time,
        tracer_time,
        "automatic source temperature replacement",
    )

    if diagnostics_dir is not None:
        diagnostics_dir = Path(diagnostics_dir)
        diagnostics_dir.mkdir(parents=True, exist_ok=True)

    inputs = _prepare_temperature_replacement_inputs(
        source_eles=source_eles,
        source_time_and_data=source_time_and_data,
        xctr=xctr,
        yctr=yctr,
        source_mapping_dir=source_mapping_dir,
        start_time=start_time,
        usgs_cache_folder=usgs_cache_folder,
        nwm_shapefile=nwm_shapefile,
        states=states,
        diagnostics_dir=diagnostics_dir,
        pooling=pooling,
        nwm_data_dir=nwm_data_dir,
        warning_log_path=warning_log_path,
        search_directions=search_directions,
        hgrid=hgrid,
    )
    replaced_count, total_usgs_records = _apply_temperature_replacements(
        source_eles=source_eles,
        temperature_data=temperature_data,
        inputs=inputs,
        pooling=pooling,
        diagnostics_dir=diagnostics_dir,
    )

    msource_data_list[0] = (tracer_time, temperature_data)

    print(
        f"[SOURCE TEMPERATURE] {pooling} replacement "
        f"completed for {replaced_count}/{len(source_eles)} eligible "
        f"source(s), {total_usgs_records} USGS-backed records total."
    )

    return msource_data_list, replaced_count


def replace_source_temperatures_with_usgs(
    base_ss,
    hgrid,
    source_mapping_dir: str | Path,
    start_time,
    usgs_cache_folder: str | Path,
    nwm_shapefile: str | Path,
    states=None,
    diagnostics_dir: str | Path | None = None,
    pooling: str = FIRST_USABLE,
    nwm_data_dir: str | Path | None = None,
    warning_log_path: str | Path | None = None,
    search_directions: tuple[str, ...] | list[str] = (
        TEMPERATURE_SEARCH_DIRECTIONS
    ),
):
    """Return a copy with eligible source temperatures replaced by USGS.

    This stage changes only the first ``msource`` tracer. Source flow, sink
    flow, element IDs, column ordering, and the input object are unchanged.
    """
    source_eles, source_values, msource_values = _copy_source_components(
        base_ss
    )
    sink_eles, sink_values = _copy_sink_components(base_ss)
    xctr, yctr = _compute_grid_centers(hgrid)

    msource_values, replaced_count = (
        _apply_automatic_temperature_replacements(
            source_eles=source_eles,
            source_time_and_data=source_values,
            msource_data_list=msource_values,
            xctr=xctr,
            yctr=yctr,
            source_mapping_dir=source_mapping_dir,
            start_time=start_time,
            usgs_cache_folder=usgs_cache_folder,
            nwm_shapefile=nwm_shapefile,
            states=list(states or STOFS3D_ATL_STATES),
            diagnostics_dir=(
                None
                if diagnostics_dir is None
                else Path(diagnostics_dir)
            ),
            pooling=pooling,
            nwm_data_dir=nwm_data_dir,
            warning_log_path=warning_log_path,
            search_directions=search_directions,
            hgrid=hgrid,
        )
    )

    corrected_ss = _build_source_sink(
        source_eles=source_eles,
        source_time_and_data=source_values,
        msource_data_list=msource_values,
        sink_eles=sink_eles,
        sink_time_and_data=sink_values,
    )
    return corrected_ss, replaced_count
