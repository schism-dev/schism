#!/usr/bin/env python3
"""
Post-process an already assembled SCHISM ``source_sink`` object for
artificial-island rivers.

This module intentionally does not call or modify:

    gen_sourcesink_nwm()
    source_nwm2usgs()
    relocate_sources2()
    set_constant_sink()

The main function operates directly on ``base_ss`` after source relocation
and the standard USGS replacement procedure:

    base_ss = apply_artificial_island_corrections(
        base_ss=base_ss,
        hgrid=hgrid,
        original_source_sink_dir=original_source_sink_dir,
        patch_info_file=config.artificial_island_source_sink_info,
        start_time=config.startdate,
        usgs_cache_folder=usgs_cache_folder,
        output_dir=Path(wdir) / "patch_artificial_island_source_sink",
    )

Supported YAML patch types
--------------------------

1. ``force_source_sink_locations``

   Create a source or sink at the main-grid element nearest the specified
   location. A nearby original forcing column is used when available.
   For source entries, USGS flow and temperature observations may optionally
   replace the original forcing.

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
       nearby match. If no original forcing lies within this radius, the
       nearest available original forcing may still be used as a fallback.
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
1. Copy the relocated source/sink forcing from ``base_ss``.
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

from pylib import schism_grid as read_schism_grid
from pylib_experimental.schism_file import source_sink
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    CFS_TO_CMS,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.usgs_series import (
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
from stofs3d_setup.ops.Source_sink.spatial_corrections import (
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
from stofs3d_setup.ops.Source_sink.source_sink_diagnostics import (
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


def apply_artificial_island_corrections(
    base_ss: source_sink,
    hgrid,
    original_source_sink_dir,
    patch_info_file: str | Path | dict,
    start_time=None,
    usgs_cache_folder=None,
    output_dir: str | Path | None = None,
) -> source_sink:
    """
    Restore only YAML-listed artificial-island source/sink forcing.

    The forcing columns are copied from ``original_source_sink_dir`` and
    appended to the already relocated ``base_ss``. Existing relocated source
    or sink columns are never moved.

    """
    if base_ss is None:
        raise ValueError("base_ss must not be None")

    if start_time is None:
        raise ValueError("start_time is required for artificial-island patching")
    usgs_cache_folder = Path(
        usgs_cache_folder or Path(output_dir or Path.cwd()) / "USGS_cache"
    )

    original_source_sink_dir = Path(original_source_sink_dir)

    if not original_source_sink_dir.is_dir():
        raise FileNotFoundError(
            "Original source/sink directory does not exist: "
            f"{original_source_sink_dir}"
        )

    original_hgrid_file = original_source_sink_dir / "hgrid.gr3"
    if not original_hgrid_file.is_file():
        raise FileNotFoundError(
            "Original source/sink grid does not exist: "
            f"{original_hgrid_file}"
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

    print(
        "[ARTIFICIAL ISLAND PATCH] reading original forcing from "
        f"{original_source_sink_dir}"
    )
    original_ss = source_sink.from_files(
        source_dir=str(original_source_sink_dir),
    )
    original_hgrid = read_schism_grid(str(original_hgrid_file))

    source_eles, source_time_and_data, msource_data_list = (
        _copy_source_components(base_ss)
    )
    sink_eles, sink_time_and_data = _copy_sink_components(base_ss)

    original_source_eles, original_source_time_and_data, original_msource = (
        _copy_source_components(original_ss)
    )
    original_sink_eles, original_sink_time_and_data = (
        _copy_sink_components(original_ss)
    )

    xctr, yctr = _compute_grid_centers(hgrid)
    original_xctr, original_yctr = _compute_grid_centers(original_hgrid)
    transformer = _make_transformer()

    print(
        "[ARTIFICIAL ISLAND PATCH] "
        f"starting relocated base: {len(source_eles)} source(s), "
        f"{len(sink_eles)} sink(s)."
    )
    print(
        "[ARTIFICIAL ISLAND PATCH] "
        f"available original forcing: {len(original_source_eles)} source(s), "
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
    # Restore only YAML-listed source/sink columns from original forcing.
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

            # Use the original nearest source as the msource template even when
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
                    start_time=start_time,
                    model_time=model_time_for_obs,
                    usgs_cache_folder=usgs_cache_folder,
                )

            if usgs_flow_series is not None:
                if not original_available:
                    raise ValueError(
                        f"[ARTIFICIAL ISLAND PATCH] {name}: USGS flow is "
                        "available, but no original source forcing exists "
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
                    "records; original NWM supplied "
                    f"{int((~use_usgs_flow).sum())} records."
                )

                if tracer_columns:
                    print(
                        f"[ARTIFICIAL ISLAND PATCH] {name}: "
                        f"USGS temperature station {temperature_station_id} "
                        f"supplied temperature msource for "
                        f"{int(use_usgs_temperature.sum())}/"
                        f"{len(use_usgs_temperature)} records; original "
                        "msource supplied "
                        f"{int((~use_usgs_temperature).sum())} records."
                    )

                source_origin_message = (
                    f"USGS flow station {flow_station_id} blended with original "
                    f"source element {original_ele} "
                    f"(distance={original_distance_m:.1f} m)"
                )

            elif original_available:
                new_time = original_time
                new_flow = template_flow
                tracer_columns = template_tracers

                if original_within_radius:
                    source_origin_message = (
                        f"original source element {original_ele} "
                        f"(distance={original_distance_m:.1f} m)"
                    )
                else:
                    source_origin_message = (
                        f"nearest original source element {original_ele} "
                        f"outside radius (distance={original_distance_m:.1f} m); "
                        "restored because USGS observation was unavailable"
                    )
                    print(
                        f"[ARTIFICIAL ISLAND PATCH] warning: {name}: "
                        f"no original source within {radius_m:.1f} m. "
                        f"Creating source at main-grid element {target_ele} "
                        f"using nearest original forcing {original_ele}."
                    )
            else:
                raise ValueError(
                    f"[ARTIFICIAL ISLAND PATCH] {name}: neither USGS flow "
                    "nor an original source forcing column is available."
                )

            (
                source_eles,
                source_time_and_data,
                msource_data_list,
            ) = _append_source_column(
                source_eles=source_eles,
                source_time_and_data=source_time_and_data,
                msource_data_list=msource_data_list,
                target_ele=target_ele,
                source_time=new_time,
                flow=new_flow,
                tracer_columns=tracer_columns,
            )

            source_time, source_data = source_time_and_data
            source_idx = len(source_eles) - 1

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
                f"[ARTIFICIAL ISLAND PATCH] {name}: created source at "
                f"main-grid element {target_ele}; target-center distance="
                f"{target_distance_m:.1f} m; forcing={source_origin_message}."
            )

        else:
            if original_sink_time_and_data is None:
                raise ValueError(
                    f"[ARTIFICIAL ISLAND PATCH] {name}: "
                    "original source/sink forcing contains no vsink"
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
                    f"[ARTIFICIAL ISLAND PATCH] {name}: original forcing "
                    "contains no sink column that can be copied."
                )

            if original_distance_m > radius_m:
                print(
                    f"[ARTIFICIAL ISLAND PATCH] warning: {name}: no "
                    f"original sink within {radius_m:.1f} m. Creating "
                    f"sink at main-grid element {target_ele} using nearest "
                    f"original sink {original_ele} at "
                    f"{original_distance_m:.1f} m."
                )

            if target_ele in sink_eles:
                raise ValueError(
                    f"[ARTIFICIAL ISLAND PATCH] {name}: target element "
                    f"{target_ele} already contains a sink. Refusing to "
                    "create duplicate sink forcing."
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
                sink_data = np.column_stack((sink_data, restored_sink))

            sink_eles.append(int(target_ele))
            sink_time_and_data = (sink_time, sink_data)
            restored_sink_count += 1

            print(
                f"[ARTIFICIAL ISLAND PATCH] {name}: restored original "
                f"sink element {original_ele} -> main-grid element "
                f"{target_ele}; original distance="
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


def zero_artificial_island_sources_after_replace_USGS_before_relocation(
    source_sink_dir: str | Path,
    hgrid_file: str | Path,
    patch_info_file: str | Path | dict,
    output_dir: str | Path | None = None,
) -> None:
    """
    Set vsource columns to zero for pre-relocation source elements located
    within 50 m of YAML-listed artificial-island locations.

    This function is intended to run after source_nwm2usgs() and before
    relocate_sources2().

    Only source_sink_dir/vsource.th is modified. Source elements, msource,
    vsink, source_sink.in, sources.json, and sinks.json remain unchanged.
    """
    search_radius_m = 500.0

    source_sink_dir = Path(source_sink_dir)
    hgrid_file = Path(hgrid_file)

    if not source_sink_dir.is_dir():
        raise FileNotFoundError(
            f"Source/sink directory does not exist: {source_sink_dir}"
        )

    if not hgrid_file.is_file():
        raise FileNotFoundError(
            f"Pre-relocation hgrid does not exist: {hgrid_file}"
        )

    patch_info = _load_patch_info(patch_info_file)

    raw_points = (
        patch_info.get(
            "remove_source_locations_in_artificial_island"
        )
        or []
    )

    if not raw_points:
        print(
            "[PRE-RELOCATION SOURCE ZERO] no artificial-island "
            "source locations were specified; vsource.th is unchanged."
        )
        return

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

    # Read the exact grid used to generate the pre-relocation source elements.
    hgrid = read_schism_grid(str(hgrid_file))

    original_ss = source_sink.from_files(
        source_dir=str(source_sink_dir),
    )

    source_eles = [
        int(ele)
        for ele in np.asarray(original_ss.source_eles).reshape(-1)
    ]

    if original_ss.vsource is None or not source_eles:
        print(
            "[PRE-RELOCATION SOURCE ZERO] no source forcing exists; "
            "nothing was changed."
        )
        return

    source_time = _time_array(original_ss.vsource)
    source_data = _data_array(original_ss.vsource).copy()

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
                "[PRE-RELOCATION SOURCE ZERO] "
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
            "[PRE-RELOCATION SOURCE ZERO] "
            f"{point['name']}: set {len(found_elements)} source(s) "
            f"to zero within {search_radius_m:.1f} m: "
            f"{sorted(found_elements)}"
        )

    # Overwrite only vsource.th.
    np.savetxt(
        source_sink_dir / "vsource.th",
        np.c_[source_time, source_data],
        fmt="%.6f",
    )

    print(
        "[PRE-RELOCATION SOURCE ZERO] completed: "
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
            / "zeroed_sources_after_replace_USGS_before_relocation.csv",
            index=False,
        )

        np.savetxt(
            output_dir / "vsource_zeroed.th",
            np.c_[source_time, source_data],
            fmt="%.6f",
        )

__all__ = [
    "apply_artificial_island_corrections",
    "zero_artificial_island_sources_after_replace_USGS_before_relocation",
]
