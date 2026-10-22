#!/usr/bin/env python3
"""
Post-process a SCHISM ``source_sink`` object for
artificial-island rivers.

This module intentionally does not call or modify:

    gen_sourcesink_nwm()
    source_nwm2usgs()
    relocate_sources2()
    set_constant_sink()

The main function uses the current forcing and a separate reference snapshot:

    base_ss = apply_artificial_island_corrections(
        base_ss=base_ss,
        hgrid=hgrid,
        reference_ss=reference_ss,
        reference_hgrid=reference_hgrid,
        patch_info_file=config.artificial_island_source_sink_info,
        start_time=config.startdate,
        usgs_cache_folder=usgs_cache_folder,
        output_dir=Path(wdir) / "patch_artificial_island_source_sink",
    )

Supported YAML patch types
--------------------------

1. ``force_source_sink_locations``

   Create a source or sink at the main-grid element nearest the specified
   location. A nearby reference forcing column is used when available.
   For source entries, USGS flow and temperature observations may optionally
   replace the reference forcing.

   Parameters
   ----------
   name : str
       Descriptive location name. When ``use_usgs_obs`` is true, the name
       must have a corresponding station entry in ``USGS_STATION_BY_NAME``.
   source_sink_type : {"source", "sink"}
       Type of forcing to create.
   x, y : float
       Longitude and latitude of the target location.
   max_search_radius_m : float
       Distance used to determine whether an original source or sink is a
       nearby match. If no reference forcing lies within this radius, the
       nearest available reference forcing may still be used as a fallback.
   allow_negative_sink : bool
       For source entries, move negative flow values to ``vsink`` at the
       same element and clip the corresponding ``vsource`` values to zero.
   use_usgs_obs : bool
       Replace available source-flow and temperature records with USGS
       observations. Original forcing is retained outside valid observation
       coverage.

   Example
   -------
   force_source_sink_locations:
     - name: Wando
       source_sink_type: source
       x: -79.6996536667
       y: 32.9802923333
       max_search_radius_m: 1500.0
       allow_negative_sink: false
       use_usgs_obs: true

2. ``large_constant_sink_artificial_island_locations``

   Add a constant sink of -1000 m3/s at the main-grid element nearest each
   listed longitude/latitude location.

   Parameters
   ----------
   name : str
       Descriptive artificial-island or river location name.
   x, y : float
       Longitude and latitude used to select the nearest main-grid element.

   Example
   -------
   large_constant_sink_artificial_island_locations:
     - name: Wando
       x: -79.6996536667
       y: 32.9802923333

Processing order
----------------
1. Copy the current source/sink forcing from ``base_ss``.
2. Add -1000 m3/s constant sinks for entries under
   ``large_constant_sink_artificial_island_locations``.
3. Restore or create entries under ``force_source_sink_locations``.
4. Return a new ``source_sink`` object. The input ``base_ss`` is not modified.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer

from pylib_experimental.schism_file import source_sink
from stofs3d_setup.ops.USGS.ogc.store import DEFAULT_USGS_CACHE_FOLDER
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    CFS_TO_CMS,
)
from stofs3d_setup.ops.USGS.ogc.series import (
    _get_usgs_forcing,
    _interpolate_usgs_with_original_fallback,
    _model_datetimes,
)
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _add_negative_source_part_to_sink,
    _append_source_column,
    _build_source_sink,
    _check_matching_time,
    _copy_sink_components,
    _copy_source_components,
    _data_array,
    _time_array,
)
from stofs3d_setup.ops.Source_sink.spatial_selection import (
    _compute_grid_centers,
    _elements_within_radius,
    _make_transformer,
    _nearest_element,
)
from stofs3d_setup.ops.Source_sink.Patch_artificial_island.island_config import (
    _as_dict,
    _load_patch_info,
    _normalize_force_points,
    _normalize_large_constant_sink_points,
)
from stofs3d_setup.ops.Source_sink.Patch_artificial_island.island_diagnostics import (
    _write_diagnostics,
)


def _add_large_constant_sinks(
    points: list[dict],
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
    source_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    sink_eles: list[int],
    sink_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    sink_value: float = -1000.0,
    diagnostics_dir: Path | None = None,
) -> tuple[
    list[int],
    tuple[np.ndarray, np.ndarray] | None,
    int,
]:
    """
    Add a constant sink at the main-grid element nearest each YAML location.

    If a sink already exists at a target element, ``sink_value`` is added to
    the existing sink column. Repeated YAML locations resolving to the same
    element therefore accumulate.
    """
    if not points:
        return sink_eles, sink_time_and_data, 0

    if sink_value >= 0.0:
        raise ValueError(
            f"sink_value must be negative; received {sink_value}"
        )

    if source_time_and_data is not None:
        model_time = np.asarray(
            source_time_and_data[0],
            dtype=float,
        ).copy()
    elif sink_time_and_data is not None:
        model_time = np.asarray(
            sink_time_and_data[0],
            dtype=float,
        ).copy()
    else:
        raise ValueError(
            "No source or sink time array is available for creating "
            "large artificial-island constant sinks"
        )

    all_element_ids = np.arange(
        1,
        len(xctr) + 1,
        dtype=int,
    )

    affected_elements: set[int] = set()
    diagnostic_rows: list[dict] = []

    for point in points:
        target_ele, distance_m = _nearest_element(
            x=point["x"],
            y=point["y"],
            candidate_element_ids=all_element_ids,
            xctr=xctr,
            yctr=yctr,
            transformer=transformer,
        )

        if target_ele is None:
            raise ValueError(
                f"[ARTIFICIAL ISLAND PATCH] {point['name']}: "
                "could not find the nearest main-grid element"
            )

        constant_values = np.full(
            len(model_time),
            float(sink_value),
            dtype=float,
        )

        if sink_time_and_data is None:
            sink_eles = [int(target_ele)]
            sink_time_and_data = (
                model_time.copy(),
                constant_values.reshape(-1, 1),
            )
            action = "created"
        else:
            sink_time, sink_data = sink_time_and_data
            _check_matching_time(
                model_time,
                sink_time,
                f"large constant sink for {point['name']}",
            )

            if target_ele in sink_eles:
                sink_idx = sink_eles.index(target_ele)
                sink_data[:, sink_idx] += constant_values
                action = "added_to_existing"
            else:
                sink_eles.append(int(target_ele))
                sink_data = np.column_stack(
                    (sink_data, constant_values)
                )
                action = "created"

            sink_time_and_data = (sink_time, sink_data)

        affected_elements.add(int(target_ele))

        diagnostic_rows.append(
            {
                "name": point["name"],
                "x": point["x"],
                "y": point["y"],
                "sink_element": int(target_ele),
                "distance_m": float(distance_m),
                "constant_sink_m3s": float(sink_value),
                "action": action,
            }
        )

        print(
            "[ARTIFICIAL ISLAND PATCH] large constant sink: "
            f"{point['name']} -> element {target_ele}, "
            f"distance={distance_m:.1f} m, "
            f"sink={sink_value:.1f} m3/s, "
            f"action={action}."
        )

    if diagnostics_dir is not None:
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            diagnostic_rows,
            columns=[
                "name",
                "x",
                "y",
                "sink_element",
                "distance_m",
                "constant_sink_m3s",
                "action",
            ],
        ).to_csv(
            diagnostics_dir
            / "large_constant_sink_artificial_island_locations.csv",
            index=False,
        )

    return (
        sink_eles,
        sink_time_and_data,
        len(affected_elements),
    )


def _set_or_append_source_column(
    source_eles,
    source_time_and_data,
    msource_data_list,
    target_ele,
    source_time,
    flow,
    tracer_columns,
):
    """Restore forcing at an existing source, or create a new source."""
    if target_ele not in source_eles:
        return _append_source_column(
            source_eles=source_eles,
            source_time_and_data=source_time_and_data,
            msource_data_list=msource_data_list,
            target_ele=target_ele,
            source_time=source_time,
            flow=flow,
            tracer_columns=tracer_columns,
        )

    source_idx = source_eles.index(target_ele)
    current_time, source_data = source_time_and_data
    _check_matching_time(current_time, source_time, f"source {target_ele}")
    flow = np.asarray(flow, dtype=float).reshape(-1)
    if len(flow) != len(current_time):
        raise ValueError(f"Source {target_ele} flow length does not match model time")
    if len(msource_data_list) != len(tracer_columns):
        raise ValueError("Reference tracer count does not match existing source")

    source_data[:, source_idx] = flow
    for tracer_idx, ((tracer_time, tracer_data), values) in enumerate(
        zip(msource_data_list, tracer_columns)
    ):
        _check_matching_time(
            tracer_time, source_time, f"source {target_ele} tracer {tracer_idx + 1}"
        )
        values = np.asarray(values, dtype=float).reshape(-1)
        if len(values) != len(tracer_time):
            raise ValueError(
                f"Source {target_ele} tracer {tracer_idx + 1} length does not match model time"
            )
        tracer_data[:, source_idx] = values

    return source_eles, source_time_and_data, msource_data_list


def apply_artificial_island_corrections(
    base_ss: source_sink,
    hgrid,
    reference_ss: source_sink,
    reference_hgrid,
    patch_info_file: str | Path | dict,
    start_time=None,
    usgs_cache_folder=None,
    output_dir: str | Path | None = None,
) -> source_sink:
    """
    Restore only YAML-listed artificial-island source/sink forcing.

    Selected forcing columns are copied from ``reference_ss`` to ``base_ss``.
    Existing target columns are replaced; otherwise columns are appended.
    Neither input is modified, and relocation is not required.

    """
    if base_ss is None:
        raise ValueError("base_ss must not be None")
    if reference_ss is None or reference_hgrid is None:
        raise ValueError("reference_ss and reference_hgrid are required")

    if start_time is None:
        raise ValueError("start_time is required for artificial-island patching")
    usgs_cache_folder = Path(
        usgs_cache_folder or DEFAULT_USGS_CACHE_FOLDER
    )

    patch_info = _load_patch_info(patch_info_file)
    unsupported_sections = [
        section
        for section in (
            "zero_source_regions",
            "exclude_source_sink_locations",
        )
        if patch_info.get(section)
    ]
    if unsupported_sections:
        raise ValueError(
            "Generic spatial corrections are not artificial-island "
            f"operations: {unsupported_sections}"
        )
    force_points = _normalize_force_points(patch_info)
    large_constant_sink_points = (
        _normalize_large_constant_sink_points(patch_info)
    )

    if (
        not force_points
        and not large_constant_sink_points
    ):
        print("[ARTIFICIAL ISLAND PATCH] no YAML entries; returning base_ss.")
        return deepcopy(base_ss)

    source_eles, source_time_and_data, msource_data_list = (
        _copy_source_components(base_ss)
    )
    sink_eles, sink_time_and_data = _copy_sink_components(base_ss)

    original_source_eles, original_source_time_and_data, original_msource = (
        _copy_source_components(reference_ss)
    )
    original_sink_eles, original_sink_time_and_data = (
        _copy_sink_components(reference_ss)
    )

    xctr, yctr = _compute_grid_centers(hgrid)
    original_xctr, original_yctr = _compute_grid_centers(reference_hgrid)
    transformer = _make_transformer()

    print(
        "[ARTIFICIAL ISLAND PATCH] "
        f"starting base: {len(source_eles)} source(s), "
        f"{len(sink_eles)} sink(s)."
    )
    print(
        "[ARTIFICIAL ISLAND PATCH] "
        f"available reference forcing: {len(original_source_eles)} source(s), "
        f"{len(original_sink_eles)} sink(s)."
    )

    restored_source_count = 0
    restored_sink_count = 0
    large_constant_sink_count = 0

    # ------------------------------------------------------------------
    # Add a -1000 m3/s constant sink at the main-grid element nearest
    # each YAML-listed artificial-island location.
    # ------------------------------------------------------------------
    if large_constant_sink_points:
        large_sink_diagnostics_dir = (
            None
            if output_dir is None
            else Path(output_dir) / "large_constant_sinks"
        )

        (
            sink_eles,
            sink_time_and_data,
            large_constant_sink_count,
        ) = _add_large_constant_sinks(
            points=large_constant_sink_points,
            xctr=xctr,
            yctr=yctr,
            transformer=transformer,
            source_time_and_data=source_time_and_data,
            sink_eles=sink_eles,
            sink_time_and_data=sink_time_and_data,
            sink_value=-1000.0,
            diagnostics_dir=large_sink_diagnostics_dir,
        )

    # ------------------------------------------------------------------
    # Restore only YAML-listed source/sink columns from reference forcing.
    # ------------------------------------------------------------------
    for point in force_points:
        name = point["name"]
        source_sink_type = point["source_sink_type"]
        radius_m = point["max_search_radius_m"]

        target_ele, target_distance_m = _nearest_element(
            x=point["x"],
            y=point["y"],
            candidate_element_ids=np.arange(1, len(xctr) + 1, dtype=int),
            xctr=xctr,
            yctr=yctr,
            transformer=transformer,
        )

        if target_ele is None:
            raise ValueError(
                f"[ARTIFICIAL ISLAND PATCH] {name}: "
                "could not determine a target element on the main grid"
            )

        if source_sink_type == "source":
            original_ele, original_distance_m = _nearest_element(
                x=point["x"],
                y=point["y"],
                candidate_element_ids=original_source_eles,
                xctr=original_xctr,
                yctr=original_yctr,
                transformer=transformer,
            )

            original_available = (
                original_ele is not None
                and original_source_time_and_data is not None
            )
            original_within_radius = (
                original_available
                and original_distance_m <= radius_m
            )

            # Use the nearest reference source as the msource template even when
            # it lies outside max_search_radius_m. The flow may still be
            # replaced entirely by the explicitly mapped USGS observation.
            if original_available:
                original_idx = original_source_eles.index(original_ele)
                original_time, original_data = original_source_time_and_data
                template_flow = original_data[:, original_idx].copy()
                template_tracers = [
                    tracer_data[:, original_idx].copy()
                    for _, tracer_data in original_msource
                ]
            else:
                original_idx = None
                original_time = (
                    source_time_and_data[0].copy()
                    if source_time_and_data is not None
                    else None
                )
                template_flow = None
                template_tracers = []

            use_usgs = bool(point.get("use_usgs_obs", False))
            usgs_flow_series = None
            usgs_temperature_series = None
            flow_station_id = None
            temperature_station_id = None

            if use_usgs:
                model_time_for_obs = (
                    source_time_and_data[0]
                    if source_time_and_data is not None
                    else original_time
                )
                if model_time_for_obs is None:
                    raise ValueError(
                        f"[ARTIFICIAL ISLAND PATCH] {name}: no model time "
                        "array is available for USGS interpolation"
                    )

                (
                    usgs_flow_series,
                    usgs_temperature_series,
                    flow_station_id,
                    temperature_station_id,
                ) = _get_usgs_forcing(
                    name=name,
                    flow_station_id=point.get("flow_station_id"),
                    temperature_station_id=point.get(
                        "temperature_station_id"
                    ),
                    start_time=start_time,
                    model_time=model_time_for_obs,
                    usgs_cache_folder=usgs_cache_folder,
                )

            if usgs_flow_series is not None:
                if not original_available:
                    raise ValueError(
                        f"[ARTIFICIAL ISLAND PATCH] {name}: USGS flow is "
                        "available, but no reference source forcing exists "
                        "to fill missing USGS periods and provide msource."
                    )

                new_time = (
                    source_time_and_data[0]
                    if source_time_and_data is not None
                    else original_time
                )
                target_datetime = _model_datetimes(
                    start_time,
                    new_time,
                )

                (
                    new_flow,
                    use_usgs_flow,
                ) = _interpolate_usgs_with_original_fallback(
                    series=usgs_flow_series,
                    target_time=target_datetime,
                    original_values=template_flow,
                    scale=CFS_TO_CMS,
                )

                tracer_columns = [
                    values.copy() for values in template_tracers
                ]

                use_usgs_temperature = np.zeros(
                    len(new_time),
                    dtype=bool,
                )
                if (
                    usgs_temperature_series is not None
                    and tracer_columns
                ):
                    (
                        tracer_columns[0],
                        use_usgs_temperature,
                    ) = _interpolate_usgs_with_original_fallback(
                        series=usgs_temperature_series,
                        target_time=target_datetime,
                        original_values=tracer_columns[0],
                        scale=1.0,
                    )

                print(
                    f"[ARTIFICIAL ISLAND PATCH] {name}: USGS flow station "
                    f"{flow_station_id} supplied vsource for "
                    f"{int(use_usgs_flow.sum())}/{len(use_usgs_flow)} "
                    "records; reference forcing supplied "
                    f"{int((~use_usgs_flow).sum())} records."
                )

                if tracer_columns:
                    print(
                        f"[ARTIFICIAL ISLAND PATCH] {name}: "
                        f"USGS temperature station {temperature_station_id} "
                        f"supplied temperature msource for "
                        f"{int(use_usgs_temperature.sum())}/"
                        f"{len(use_usgs_temperature)} records; reference "
                        "msource supplied "
                        f"{int((~use_usgs_temperature).sum())} records."
                    )

                source_origin_message = (
                    f"USGS flow station {flow_station_id} blended with reference "
                    f"source element {original_ele} "
                    f"(distance={original_distance_m:.1f} m)"
                )

            elif original_available:
                new_time = original_time
                new_flow = template_flow
                tracer_columns = template_tracers

                if original_within_radius:
                    source_origin_message = (
                        f"reference source element {original_ele} "
                        f"(distance={original_distance_m:.1f} m)"
                    )
                else:
                    source_origin_message = (
                        f"nearest reference source element {original_ele} "
                        f"outside radius (distance={original_distance_m:.1f} m); "
                        "restored because USGS observation was unavailable"
                    )
                    print(
                        f"[ARTIFICIAL ISLAND PATCH] warning: {name}: "
                        f"no reference source within {radius_m:.1f} m. "
                        f"Creating source at main-grid element {target_ele} "
                        f"using nearest reference forcing {original_ele}."
                    )
            else:
                raise ValueError(
                    f"[ARTIFICIAL ISLAND PATCH] {name}: neither USGS flow "
                    "nor a reference source forcing column is available."
                )

            (
                source_eles,
                source_time_and_data,
                msource_data_list,
            ) = _set_or_append_source_column(
                source_eles=source_eles,
                source_time_and_data=source_time_and_data,
                msource_data_list=msource_data_list,
                target_ele=target_ele,
                source_time=new_time,
                flow=new_flow,
                tracer_columns=tracer_columns,
            )

            source_time, source_data = source_time_and_data
            source_idx = source_eles.index(target_ele)

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
                        f"[ARTIFICIAL ISLAND PATCH] {name}: split "
                        f"{n_negative} negative record(s) into vsink at "
                        f"element {target_ele}; minimum="
                        f"{min_negative:.6f} m3/s."
                    )
                else:
                    print(
                        f"[ARTIFICIAL ISLAND PATCH] {name}: "
                        "allow_negative_sink enabled; no negative records."
                    )

            restored_source_count += 1
            print(
                f"[ARTIFICIAL ISLAND PATCH] {name}: set source at "
                f"main-grid element {target_ele}; target-center distance="
                f"{target_distance_m:.1f} m; forcing={source_origin_message}."
            )

        else:
            if original_sink_time_and_data is None:
                raise ValueError(
                    f"[ARTIFICIAL ISLAND PATCH] {name}: "
                    "reference source/sink forcing contains no vsink"
                )

            original_ele, original_distance_m = _nearest_element(
                x=point["x"],
                y=point["y"],
                candidate_element_ids=original_sink_eles,
                xctr=original_xctr,
                yctr=original_yctr,
                transformer=transformer,
            )

            if original_ele is None:
                raise ValueError(
                    f"[ARTIFICIAL ISLAND PATCH] {name}: reference forcing "
                    "contains no sink column that can be copied."
                )

            if original_distance_m > radius_m:
                print(
                    f"[ARTIFICIAL ISLAND PATCH] warning: {name}: no "
                    f"reference sink within {radius_m:.1f} m. Creating "
                    f"sink at main-grid element {target_ele} using nearest "
                    f"reference sink {original_ele} at "
                    f"{original_distance_m:.1f} m."
                )

            original_idx = original_sink_eles.index(original_ele)
            original_time, original_data = original_sink_time_and_data
            restored_sink = original_data[:, original_idx].copy()

            if sink_time_and_data is None:
                sink_time = original_time.copy()
                sink_data = restored_sink.reshape(-1, 1)
            else:
                sink_time, sink_data = sink_time_and_data
                _check_matching_time(
                    sink_time,
                    original_time,
                    f"sink {name}",
                )
                if target_ele in sink_eles:
                    sink_data[:, sink_eles.index(target_ele)] = restored_sink
                else:
                    sink_data = np.column_stack((sink_data, restored_sink))

            if target_ele not in sink_eles:
                sink_eles.append(int(target_ele))
            sink_time_and_data = (sink_time, sink_data)
            restored_sink_count += 1

            print(
                f"[ARTIFICIAL ISLAND PATCH] {name}: restored reference "
                f"sink element {original_ele} -> main-grid element "
                f"{target_ele}; reference distance="
                f"{original_distance_m:.1f} m, target-center distance="
                f"{target_distance_m:.1f} m."
            )

    if len(source_eles) != len(set(source_eles)):
        raise ValueError(
            "Artificial-island patch produced duplicate source elements"
        )
    if len(sink_eles) != len(set(sink_eles)):
        raise ValueError(
            "Artificial-island patch produced duplicate sink elements"
        )

    patched_ss = _build_source_sink(
        source_eles=source_eles,
        source_time_and_data=source_time_and_data,
        msource_data_list=msource_data_list,
        sink_eles=sink_eles,
        sink_time_and_data=sink_time_and_data,
    )

    print(
        "[ARTIFICIAL ISLAND PATCH] added large constant sinks to "
        f"{large_constant_sink_count} unique artificial-island element(s)."
    )
    print(
        "[ARTIFICIAL ISLAND PATCH] restored "
        f"{restored_source_count} YAML source(s) and "
        f"{restored_sink_count} YAML sink(s)."
    )
    print(
        "[ARTIFICIAL ISLAND PATCH] finished combined base: "
        f"{len(source_eles)} source(s), {len(sink_eles)} sink(s)."
    )

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        patched_ss.writer(str(output_dir))
        _write_diagnostics(output_dir, patched_ss, hgrid)

    return patched_ss


def suppress_artificial_island_sources(
    base_ss: source_sink,
    hgrid,
    patch_info_file: str | Path | dict,
    output_dir: str | Path | None = None,
) -> source_sink:
    """
    Return forcing with source flow zeroed near YAML-listed island locations.

    The input object is not modified. Source elements, tracers, and sinks
    remain unchanged in the returned object.
    """
    search_radius_m = 500.0
    if base_ss is None:
        raise ValueError("base_ss must not be None")

    patch_info = _load_patch_info(patch_info_file)

    raw_points = (
        patch_info.get(
            "remove_source_locations_in_artificial_island"
        )
        or []
    )

    if not raw_points:
        print(
            "[ISLAND SOURCE SUPPRESSION] no source locations specified."
        )
        return base_ss

    points: list[dict] = []

    for raw_point in raw_points:
        point = _as_dict(raw_point)

        if "x" not in point or "y" not in point:
            raise ValueError(
                "Each remove_source_locations_in_artificial_island "
                f"entry requires x and y: {point}"
            )

        points.append(
            {
                "name": str(point.get("name", "unnamed")),
                "x": float(point["x"]),
                "y": float(point["y"]),
            }
        )

    source_eles = [
        int(ele)
        for ele in np.asarray(base_ss.source_eles).reshape(-1)
    ]

    if base_ss.vsource is None or not source_eles:
        print(
            "[ISLAND SOURCE SUPPRESSION] no source forcing exists; "
            "nothing was changed."
        )
        return base_ss

    source_data = _data_array(base_ss.vsource).copy()

    if source_data.shape[1] != len(source_eles):
        raise ValueError(
            f"vsource has {source_data.shape[1]} data columns, but "
            f"source_sink.in contains {len(source_eles)} source elements"
        )

    xctr, yctr = _compute_grid_centers(hgrid)
    transformer = _make_transformer()

    source_index = {
        source_ele: source_idx
        for source_idx, source_ele in enumerate(source_eles)
    }

    zeroed_elements: set[int] = set()
    diagnostic_rows: list[dict] = []

    for point in points:
        found_elements = _elements_within_radius(
            x=point["x"],
            y=point["y"],
            radius_m=search_radius_m,
            candidate_element_ids=source_eles,
            xctr=xctr,
            yctr=yctr,
            transformer=transformer,
        )

        if not found_elements:
            nearest_ele, nearest_distance_m = _nearest_element(
                x=point["x"],
                y=point["y"],
                candidate_element_ids=source_eles,
                xctr=xctr,
                yctr=yctr,
                transformer=transformer,
            )

            nearest_text = (
                "none"
                if nearest_ele is None
                else (
                    f"element {nearest_ele} at "
                    f"{nearest_distance_m:.1f} m"
                )
            )

            print(
                "[ISLAND SOURCE SUPPRESSION] "
                f"{point['name']}: no source found within "
                f"{search_radius_m:.1f} m; nearest={nearest_text}."
            )

            diagnostic_rows.append(
                {
                    "name": point["name"],
                    "x": point["x"],
                    "y": point["y"],
                    "source_element": "",
                    "distance_m": (
                        nearest_distance_m
                        if nearest_ele is not None
                        else np.nan
                    ),
                    "mean_vsource_before": np.nan,
                    "mean_vsource_after": np.nan,
                    "action": "no_source_within_radius",
                }
            )
            continue

        point_x_m, point_y_m = transformer.transform(
            point["x"],
            point["y"],
        )

        for source_ele in found_elements:
            source_ele = int(source_ele)
            source_idx = source_index[source_ele]

            mean_before = float(
                np.mean(source_data[:, source_idx])
            )

            source_data[:, source_idx] = 0.0
            zeroed_elements.add(source_ele)

            source_x_m, source_y_m = transformer.transform(
                float(xctr[source_ele - 1]),
                float(yctr[source_ele - 1]),
            )

            distance_m = float(
                np.hypot(
                    source_x_m - point_x_m,
                    source_y_m - point_y_m,
                )
            )

            diagnostic_rows.append(
                {
                    "name": point["name"],
                    "x": point["x"],
                    "y": point["y"],
                    "source_element": source_ele,
                    "distance_m": distance_m,
                    "mean_vsource_before": mean_before,
                    "mean_vsource_after": 0.0,
                    "action": "vsource_zeroed",
                }
            )

        print(
            "[ISLAND SOURCE SUPPRESSION] "
            f"{point['name']}: set {len(found_elements)} source(s) "
            f"to zero within {search_radius_m:.1f} m: "
            f"{sorted(found_elements)}"
        )

    suppressed_ss = deepcopy(base_ss)
    suppressed_ss.vsource.df.iloc[:, :] = source_data

    print(
        "[ISLAND SOURCE SUPPRESSION] completed: "
        f"{len(zeroed_elements)} unique source element(s) set to zero."
    )

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        pd.DataFrame(
            diagnostic_rows,
            columns=[
                "name",
                "x",
                "y",
                "source_element",
                "distance_m",
                "mean_vsource_before",
                "mean_vsource_after",
                "action",
            ],
        ).to_csv(
            output_dir
            / "suppressed_artificial_island_sources.csv",
            index=False,
        )

    return suppressed_ss

__all__ = [
    "apply_artificial_island_corrections",
    "suppress_artificial_island_sources",
]
