"""Explicit flow and temperature overrides for existing SCHISM sources."""

import numpy as np
from pyproj import Transformer

from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    CFS_TO_CMS,
    USGS_FLOW_STATION_BY_NAME,
    USGS_TEMPERATURE_STATION_BY_NAME,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.usgs_series import (
    _get_usgs_flow,
    _get_usgs_temperature,
    _interpolate_usgs_with_original_fallback,
    _model_datetimes,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.source_temperature import (
    DISCHARGE_WEIGHTED,
    FIRST_USABLE,
    _blend_complete_temperature_or_preserve,
    _mixed_temperature_columns,
    _validate_temperature_pooling,
)
from stofs3d_setup.ops.Source_sink.correction_config import load_source_sink_corrections
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.override_config import (
    source_override_points,
)
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _add_negative_source_part_to_sink,
    _build_source_sink,
    _check_matching_time,
    _copy_sink_components,
    _copy_source_components,
)
from stofs3d_setup.ops.Source_sink.spatial_selection import (
    _compute_grid_centers,
    _make_transformer,
    _nearest_element,
)


def _replace_existing_relocated_source(
    point: dict,
    source_eles: list[int],
    source_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    msource_data_list: list[tuple[np.ndarray, np.ndarray]],
    sink_eles: list[int],
    sink_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
    start_time,
    usgs_cache_folder,
    temperature_pooling: str = FIRST_USABLE,
) -> tuple[
    tuple[np.ndarray, np.ndarray] | None,
    list[tuple[np.ndarray, np.ndarray]],
    list[int],
    tuple[np.ndarray, np.ndarray] | None,
]:
    """
    Replace flow and/or temperature in one existing relocated source column.

    The source element and source-column count are unchanged. Flow and
    temperature may come from different explicitly configured USGS stations.
    In ``first_usable`` mode, unsupported periods retain the existing
    relocated forcing. In ``discharge_weighted`` mode, a temperature
    override is atomic: it is applied only if every model time receives a
    finite numeric temperature.
    """
    name = point["name"]
    radius_m = point["max_search_radius_m"]
    temperature_pooling = _validate_temperature_pooling(temperature_pooling)

    if source_time_and_data is None or not source_eles:
        raise ValueError(
            f"[SOURCE OVERRIDE] {name}: relocated base has no source"
        )

    target_ele, distance_m = _nearest_element(
        x=point["x"],
        y=point["y"],
        candidate_element_ids=source_eles,
        xctr=xctr,
        yctr=yctr,
        transformer=transformer,
    )

    if target_ele is None or distance_m > radius_m:
        nearest_text = (
            "none"
            if target_ele is None
            else f"element {target_ele} at {distance_m:.1f} m"
        )
        raise ValueError(
            f"[SOURCE OVERRIDE] {name}: no relocated source found "
            f"within {radius_m:.1f} m; nearest={nearest_text}. "
            "Replace-only entries never create a new source."
        )

    flow_station_id = USGS_FLOW_STATION_BY_NAME.get(name)
    temperature_station_id = USGS_TEMPERATURE_STATION_BY_NAME.get(name)

    if point["replace_flow"] and not flow_station_id:
        raise ValueError(
            f"[SOURCE OVERRIDE] {name}: "
            "no USGS flow station is configured"
        )

    if point["replace_temperature"] and not temperature_station_id:
        print(
            f"[SOURCE OVERRIDE] warning: {name}: "
            "no USGS temperature station is configured; "
            "temperature replacement will be skipped."
        )

    source_time, source_data = source_time_and_data
    source_idx = source_eles.index(target_ele)
    target_datetime = _model_datetimes(start_time, source_time)

    flow_series = None
    temperature_series = None
    download_policy = point.get("usgs_download")
    if point["replace_flow"]:
        flow_series, flow_station_id = _get_usgs_flow(
            name=name,
            start_time=start_time,
            model_time=source_time,
            usgs_cache_folder=usgs_cache_folder,
            download_policy=download_policy,
        )
    if point["replace_temperature"]:
        temperature_series, temperature_station_id = _get_usgs_temperature(
            name=name,
            start_time=start_time,
            model_time=source_time,
            usgs_cache_folder=usgs_cache_folder,
            download_policy=download_policy,
        )

    if point["replace_flow"]:
        if flow_series is None:
            print(
                f"[SOURCE OVERRIDE] warning: {name}: "
                f"USGS flow station {flow_station_id} returned no flow; "
                "relocated vsource is unchanged."
            )
        else:
            original_flow = source_data[:, source_idx].copy()
            replaced_flow, use_usgs_flow = (
                _interpolate_usgs_with_original_fallback(
                    series=flow_series,
                    target_time=target_datetime,
                    original_values=original_flow,
                    scale=CFS_TO_CMS,
                )
            )
            source_data[:, source_idx] = replaced_flow

            print(
                f"[SOURCE OVERRIDE] {name}: replaced relocated "
                f"vsource at element {target_ele} using USGS flow station "
                f"{flow_station_id} for "
                f"{int(use_usgs_flow.sum())}/{len(use_usgs_flow)} records; "
                f"existing relocated forcing retained for "
                f"{int((~use_usgs_flow).sum())} records."
            )

            if not np.any(use_usgs_flow):
                print(
                    f"[SOURCE OVERRIDE] warning: {name}: "
                    f"USGS flow station {flow_station_id} was downloaded, "
                    "but zero model-time records were replaced. "
                    f"Model period={target_datetime[0]} to {target_datetime[-1]}; "
                    f"USGS period={flow_series.index.min()} to "
                    f"{flow_series.index.max()}."
                )

    if point["replace_temperature"]:
        if not msource_data_list:
            print(
                f"[SOURCE OVERRIDE] warning: {name}: "
                "base_ss has no msource tracer; temperature replacement "
                "was skipped."
            )
        elif temperature_station_id is None:
            print(
                f"[SOURCE OVERRIDE] warning: {name}: "
                "no USGS temperature station configured; relocated "
                "temperature msource is unchanged."
            )
        elif temperature_series is None:
            print(
                f"[SOURCE OVERRIDE] warning: {name}: "
                f"USGS temperature station {temperature_station_id} "
                "returned no temperature; relocated temperature msource "
                "is unchanged."
            )
        else:
            tracer_time, tracer_data = msource_data_list[0]
            _check_matching_time(
                source_time,
                tracer_time,
                f"relocated temperature replacement {name}",
            )

            original_temperature = tracer_data[:, source_idx].copy()
            if temperature_pooling == DISCHARGE_WEIGHTED:
                (
                    replaced_temperature,
                    use_usgs_temperature,
                    temperature_applied,
                ) = _blend_complete_temperature_or_preserve(
                    series=temperature_series,
                    target_time=target_datetime,
                    original_values=original_temperature,
                )
            else:
                replaced_temperature, use_usgs_temperature = (
                    _interpolate_usgs_with_original_fallback(
                        series=temperature_series,
                        target_time=target_datetime,
                        original_values=original_temperature,
                        scale=1.0,
                    )
                )
                temperature_applied = bool(
                    np.any(use_usgs_temperature)
                )

            tracer_data[:, source_idx] = replaced_temperature
            msource_data_list[0] = (tracer_time, tracer_data)

            if temperature_applied:
                print(
                    f"[SOURCE OVERRIDE] {name}: replaced relocated "
                    f"temperature msource tracer 1 at element {target_ele} "
                    "using USGS temperature station "
                    f"{temperature_station_id} for "
                    f"{int(use_usgs_temperature.sum())}/"
                    f"{len(use_usgs_temperature)} records; existing "
                    "temperature retained for "
                    f"{int((~use_usgs_temperature).sum())} records."
                )
            else:
                print(
                    f"[SOURCE OVERRIDE] warning: {name}: retained the "
                    "complete original temperature column because station "
                    f"{temperature_station_id} could not provide a complete "
                    "numeric replacement."
                )

    if point["allow_negative_sink"]:
        (
            sink_eles,
            sink_time_and_data,
            n_negative,
            min_negative,
        ) = _add_negative_source_part_to_sink(
            target_ele=target_ele,
            source_column_idx=source_idx,
            source_data=source_data,
            source_time=source_time,
            sink_eles=sink_eles,
            sink_time_and_data=sink_time_and_data,
        )

        if n_negative > 0:
            print(
                f"[SOURCE OVERRIDE] {name}: moved "
                f"{n_negative} negative relocated-source record(s) to vsink "
                f"at element {target_ele}; minimum={min_negative:.6f} m3/s."
            )
        else:
            print(
                f"[SOURCE OVERRIDE] {name}: "
                "allow_negative_sink enabled; no negative records."
            )

    print(
        f"[SOURCE OVERRIDE] {name}: replace-only operation matched "
        f"relocated source element {target_ele} at {distance_m:.1f} m."
    )

    return (
        (source_time, source_data),
        msource_data_list,
        sink_eles,
        sink_time_and_data,
    )


def _apply_source_overrides(
    base_ss,
    hgrid,
    points: list[dict],
    start_time,
    usgs_cache_folder,
    temperature_pooling: str = FIRST_USABLE,
):
    source_eles, source_values, msource_values = _copy_source_components(
        base_ss
    )
    sink_eles, sink_values = _copy_sink_components(base_ss)
    xctr, yctr = _compute_grid_centers(hgrid)
    transformer = _make_transformer()

    for point in points:
        source_values, msource_values, sink_eles, sink_values = (
            _replace_existing_relocated_source(
                point=point,
                source_eles=source_eles,
                source_time_and_data=source_values,
                msource_data_list=msource_values,
                sink_eles=sink_eles,
                sink_time_and_data=sink_values,
                xctr=xctr,
                yctr=yctr,
                transformer=transformer,
                start_time=start_time,
                usgs_cache_folder=usgs_cache_folder,
                temperature_pooling=temperature_pooling,
            )
        )

    if temperature_pooling == DISCHARGE_WEIGHTED and msource_values:
        mixed_columns = np.flatnonzero(
            _mixed_temperature_columns(msource_values[0][1])
        )
        if len(mixed_columns) > 0:
            mixed_elements = [source_eles[idx] for idx in mixed_columns]
            raise ValueError(
                "discharge_weighted selected-source temperature overrides "
                "left columns that are not uniformly all-ambient or "
                f"all-finite-numeric for source elements {mixed_elements}"
            )

    corrected_ss = _build_source_sink(
        source_eles=source_eles,
        source_time_and_data=source_values,
        msource_data_list=msource_values,
        sink_eles=sink_eles,
        sink_time_and_data=sink_values,
    )
    return corrected_ss, len(points)


def apply_source_flow_overrides(
    base_ss,
    hgrid,
    points: list[dict],
    start_time,
    usgs_cache_folder,
):
    """Replace configured source flows without changing source temperature."""
    flow_points = [
        {**point, "replace_temperature": False}
        for point in points
        if point.get("replace_flow", False)
    ]
    return _apply_source_overrides(
        base_ss=base_ss,
        hgrid=hgrid,
        points=flow_points,
        start_time=start_time,
        usgs_cache_folder=usgs_cache_folder,
        temperature_pooling=FIRST_USABLE,
    )


def apply_source_temperature_overrides(
    base_ss,
    hgrid,
    points: list[dict],
    start_time,
    usgs_cache_folder,
    temperature_pooling: str = FIRST_USABLE,
):
    """Replace configured source temperatures without changing source flow."""
    temperature_points = [
        {
            **point,
            "replace_flow": False,
            "allow_negative_sink": False,
        }
        for point in points
        if point.get("replace_temperature", False)
    ]
    return _apply_source_overrides(
        base_ss=base_ss,
        hgrid=hgrid,
        points=temperature_points,
        start_time=start_time,
        usgs_cache_folder=usgs_cache_folder,
        temperature_pooling=temperature_pooling,
    )


def apply_source_overrides(
    base_ss,
    hgrid,
    points: list[dict],
    start_time,
    usgs_cache_folder,
    temperature_pooling: str = FIRST_USABLE,
):
    """Apply configured USGS flow/temperature overrides to existing sources."""
    return _apply_source_overrides(
        base_ss=base_ss,
        hgrid=hgrid,
        points=points,
        start_time=start_time,
        usgs_cache_folder=usgs_cache_folder,
        temperature_pooling=temperature_pooling,
    )


def load_source_override_points(correction_info) -> list[dict]:
    """Load and normalize selected-source corrections for this stage."""
    corrections = load_source_sink_corrections(correction_info)
    return source_override_points(corrections)
