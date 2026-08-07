#!/usr/bin/env python3
"""
Post-process an already assembled SCHISM ``source_sink`` object for
artificial-island rivers and selected relocated sources.

This module intentionally does not call or modify:

    gen_sourcesink_nwm()
    source_nwm2usgs()
    relocate_sources2()
    set_constant_sink()

The main function operates directly on ``base_ss`` after source relocation
and the standard USGS replacement procedure:

    base_ss = patch_artificial_island_source_sink(
        base_ss=base_ss,
        hgrid=hgrid,
        original_source_sink_dir=original_source_sink_dir,
        relocated_source_sink_dir=relocated_source_sink_dir,
        patch_info_file=config.artificial_island_source_sink_info,
        start_time=config.startdate,
        usgs_cache_folder=usgs_cache_folder,
        nwm_shapefile=nwm_shapefile,
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

2. ``replace_only_source_locations``

   Replace flow and/or temperature values in an existing relocated source
   column. This operation does not create, remove, or move a source element.

   Parameters
   ----------
   name : str
       Descriptive location name associated with a station in
       ``USGS_STATION_BY_NAME``.
   x, y : float
       Longitude and latitude used to locate the nearest relocated source.
   max_search_radius_m : float
       Maximum permitted distance between the specified location and the
       matched relocated source element.
   replace_flow : bool
       Replace available ``vsource`` records with USGS streamflow.
   replace_temperature : bool
       Replace available temperature values in msource tracer 1.
   allow_negative_sink : bool
       Move negative replaced-flow values to ``vsink`` at the same element
       and clip the corresponding ``vsource`` values to zero.

   Example
   -------
   replace_only_source_locations:
     - name: Delaware
       x: -74.9432325
       y: 40.3422775
       max_search_radius_m: 500.0
       replace_flow: true
       replace_temperature: true
       allow_negative_sink: false

3. ``large_constant_sink_artificial_island_locations``

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

4. ``exclude_source_sink_locations``

   Remove existing source and/or sink columns whose main-grid element
   centers lie within the specified radius.

   Parameters
   ----------
   name : str
       Descriptive exclusion name.
   x, y : float
       Longitude and latitude of the exclusion center.
   radius_m : float
       Search radius around the specified coordinate.
   remove : list[str]
       Forcing types to remove. Entries may include ``source``, ``sink``,
       or both.

   Example
   -------
   exclude_source_sink_locations:
     - name: Example open boundary
       x: -79.0
       y: 33.0
       radius_m: 1000.0
       remove:
         - source
         - sink

Processing order
----------------
1. Copy the relocated source/sink forcing from ``base_ss``.
2. Optionally replace temperatures for all relocated sources using upstream
   NWM-to-USGS associations.
3. Apply explicit ``replace_only_source_locations`` replacements.
4. Add -1000 m3/s constant sinks for entries under
   ``large_constant_sink_artificial_island_locations``.
5. Restore or create entries under ``force_source_sink_locations``.
6. Remove entries under ``exclude_source_sink_locations``.
7. Return a new ``source_sink`` object. The input ``base_ss`` is not modified.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import json
from typing import Any

import numpy as np
import pandas as pd
import yaml
from pyproj import Transformer
from scipy.spatial import cKDTree

from pylib import schism_grid as read_schism_grid
from pylib_experimental.schism_file import source_sink, TimeHistory
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
from stofs3d_setup.utils.utils import STOFS3D_ATL_STATES

# USGS station mapping used when YAML has ``use_usgs_obs: true``.
# Wando is left as None until a station is selected.
USGS_STATION_BY_NAME = {
    # Forced/restored artificial-island sources
    "Wando": None,
    "Turkey": "02172035",
    "Buffalo Bluff": "02244040",
    "Dunns Creek": "02244440",

    # Existing relocated sources: replace values only
    "Delaware": "01463500",
    "Hudson River": "01358000",
}

USGS_FLOW_PARAMETER_ID = "00060"
USGS_TEMPERATURE_PARAMETER_ID = "00010"
CFS_TO_CMS = 0.028316846592
MAX_USGS_INTERPOLATION_GAP = pd.Timedelta("6 hours")


def _as_dict(value: Any) -> dict:
    """Convert a Pydantic model or mapping-like value to a plain dictionary."""
    if hasattr(value, "model_dump"):
        return value.model_dump()
    return dict(value)


def _load_patch_info(patch_info_file: str | Path | dict) -> dict:
    """Read patch settings from YAML, or accept an already loaded dictionary."""
    if isinstance(patch_info_file, dict):
        return deepcopy(patch_info_file)

    patch_info_file = Path(patch_info_file)
    if not patch_info_file.exists():
        raise FileNotFoundError(
            f"Artificial-island source/sink YAML does not exist: {patch_info_file}"
        )

    with patch_info_file.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}

    if not isinstance(data, dict):
        raise ValueError(
            f"Artificial-island YAML must contain a mapping at the top level: "
            f"{patch_info_file}"
        )

    return data


def _normalize_force_points(patch_info: dict) -> list[dict]:
    """Normalize forced source/sink locations from supported YAML keys."""
    raw_points = (
        patch_info.get("force_source_sink_locations")
        or patch_info.get("artificial_island_source_sink_locations")
        or patch_info.get("locations")
        or []
    )

    points = []
    for raw in raw_points:
        p = _as_dict(raw)

        if "x" not in p or "y" not in p:
            raise ValueError(
                f"Forced artificial-island entry requires x and y: {p}"
            )

        source_sink_type = p.get(
            "source_sink_type",
            p.get("type", p.get("kind", "source")),
        )
        source_sink_type = str(source_sink_type).lower()

        if source_sink_type not in {"source", "sink"}:
            raise ValueError(
                f"Unsupported source_sink_type={source_sink_type!r} "
                f"for entry {p.get('name', '')!r}"
            )

        points.append(
            {
                **p,
                "name": str(p.get("name", "unnamed")),
                "x": float(p["x"]),
                "y": float(p["y"]),
                "source_sink_type": source_sink_type,
                "max_search_radius_m": float(
                    p.get("max_search_radius_m", p.get("radius_m", 1500.0))
                ),
                "allow_negative_sink": bool(
                    p.get("allow_negative_sink", False)
                ),
                "use_usgs_obs": bool(p.get("use_usgs_obs", False)),
            }
        )

    return points



def _normalize_replace_relocated_points(patch_info: dict) -> list[dict]:
    """
    Normalize existing relocated-source replacement locations.

    These entries search only ``base_ss.source_eles`` and never create,
    remove, or move a source element.
    """
    raw_points = (
        patch_info.get("replace_relocated_source_locations")
        or patch_info.get("replace_only_source_locations")
        or []
    )

    points = []
    for raw in raw_points:
        p = _as_dict(raw)

        if "x" not in p or "y" not in p:
            raise ValueError(
                f"Relocated-source replacement entry requires x and y: {p}"
            )

        name = str(p.get("name", "unnamed"))
        if name not in USGS_STATION_BY_NAME:
            raise ValueError(
                f"{name!r} is listed for relocated-source replacement, "
                "but it is missing from USGS_STATION_BY_NAME"
            )

        points.append(
            {
                **p,
                "name": name,
                "x": float(p["x"]),
                "y": float(p["y"]),
                "max_search_radius_m": float(
                    p.get("max_search_radius_m", p.get("radius_m", 500.0))
                ),
                "replace_flow": bool(p.get("replace_flow", True)),
                "replace_temperature": bool(
                    p.get("replace_temperature", True)
                ),
                "allow_negative_sink": bool(
                    p.get("allow_negative_sink", False)
                ),
            }
        )

    return points



def _normalize_large_constant_sink_points(
    patch_info: dict,
) -> list[dict]:
    """Normalize artificial-island locations for large constant sinks."""
    raw_points = (
        patch_info.get(
            "large_constant_sink_artificial_island_locations"
        )
        or []
    )

    points = []
    for raw in raw_points:
        point = _as_dict(raw)

        if "x" not in point or "y" not in point:
            raise ValueError(
                "Large artificial-island constant-sink entry "
                f"requires x and y: {point}"
            )

        points.append(
            {
                "name": str(point.get("name", "unnamed")),
                "x": float(point["x"]),
                "y": float(point["y"]),
            }
        )

    return points


def _normalize_exclude_points(patch_info: dict) -> list[dict]:
    """Normalize source/sink exclusion locations."""
    raw_points = patch_info.get("exclude_source_sink_locations") or []

    points = []
    for raw in raw_points:
        p = _as_dict(raw)

        if "x" not in p or "y" not in p:
            raise ValueError(
                f"Artificial-island exclusion entry requires x and y: {p}"
            )

        remove = p.get("remove", ["source", "sink"])
        if isinstance(remove, str):
            remove = [remove]
        remove = [str(v).lower() for v in remove]

        invalid = sorted(set(remove) - {"source", "sink"})
        if invalid:
            raise ValueError(
                f"Unsupported exclusion types {invalid} for "
                f"entry {p.get('name', '')!r}"
            )

        points.append(
            {
                **p,
                "name": str(p.get("name", "unnamed")),
                "x": float(p["x"]),
                "y": float(p["y"]),
                "radius_m": float(
                    p.get("radius_m", p.get("max_search_radius_m", 500.0))
                ),
                "remove": remove,
            }
        )

    return points


def _compute_grid_centers(hgrid) -> tuple[np.ndarray, np.ndarray]:
    """
    Return element-center coordinates without leaving new xctr/yctr
    attributes on the caller's hgrid object.
    """
    had_xctr = hasattr(hgrid, "xctr")
    had_yctr = hasattr(hgrid, "yctr")

    old_xctr = (
        np.asarray(hgrid.xctr).copy()
        if had_xctr
        else None
    )
    old_yctr = (
        np.asarray(hgrid.yctr).copy()
        if had_yctr
        else None
    )

    try:
        hgrid.compute_ctr()

        xctr = np.asarray(hgrid.xctr, dtype=float).copy()
        yctr = np.asarray(hgrid.yctr, dtype=float).copy()

    finally:
        if had_xctr:
            hgrid.xctr = old_xctr
        elif hasattr(hgrid, "xctr"):
            delattr(hgrid, "xctr")

        if had_yctr:
            hgrid.yctr = old_yctr
        elif hasattr(hgrid, "yctr"):
            delattr(hgrid, "yctr")

    return xctr, yctr

def _make_transformer() -> Transformer:
    """
    Build a lon/lat-to-meter transformer.

    EPSG:3857 is used only for local nearest-neighbor screening. For the
    O(1 km) search radii used here, it is adequate and consistent with the
    earlier user-defined patch utilities.
    """
    return Transformer.from_crs(
        "EPSG:4326",
        "EPSG:3857",
        always_xy=True,
    )


def _element_center_tree(
    element_ids: list[int] | np.ndarray,
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
) -> tuple[np.ndarray, cKDTree | None]:
    """Build a KDTree for selected SCHISM element centers."""
    element_ids = np.asarray(element_ids, dtype=int).reshape(-1)

    if len(element_ids) == 0:
        return element_ids, None

    xp, yp = transformer.transform(
        xctr[element_ids - 1],
        yctr[element_ids - 1],
    )

    return element_ids, cKDTree(np.c_[xp, yp])


def _nearest_element(
    x: float,
    y: float,
    candidate_element_ids: list[int] | np.ndarray,
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
) -> tuple[int | None, float]:
    """Find the nearest candidate element center and return distance in meters."""
    element_ids, tree = _element_center_tree(
        candidate_element_ids,
        xctr,
        yctr,
        transformer,
    )

    if tree is None:
        return None, np.inf

    xq, yq = transformer.transform(x, y)
    distance_m, idx = tree.query([xq, yq])

    return int(element_ids[int(idx)]), float(distance_m)


def _elements_within_radius(
    x: float,
    y: float,
    radius_m: float,
    candidate_element_ids: list[int] | np.ndarray,
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
) -> list[int]:
    """Return candidate element IDs whose centers are within radius_m."""
    element_ids, tree = _element_center_tree(
        candidate_element_ids,
        xctr,
        yctr,
        transformer,
    )

    if tree is None:
        return []

    xq, yq = transformer.transform(x, y)
    idxs = tree.query_ball_point([xq, yq], r=float(radius_m))

    return [int(element_ids[int(idx)]) for idx in idxs]


def _time_array(th: TimeHistory) -> np.ndarray:
    """Return TimeHistory time values as a one-dimensional array."""
    return np.asarray(th.time).reshape(-1)


def _data_array(th: TimeHistory) -> np.ndarray:
    """Return TimeHistory values as a two-dimensional float array."""
    data = np.asarray(th.data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)
    return data


def _build_timehistory(
    time: np.ndarray,
    data: np.ndarray,
    element_ids: list[int],
) -> TimeHistory | None:
    """Create a TimeHistory object, or return None when no elements remain."""
    element_ids = [int(ele) for ele in element_ids]

    if len(element_ids) == 0:
        return None

    data = np.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)

    if data.shape[1] != len(element_ids):
        raise ValueError(
            "TimeHistory data/element mismatch: "
            f"{data.shape[1]} columns for {len(element_ids)} elements"
        )

    return TimeHistory(
        data_array=np.c_[np.asarray(time).reshape(-1), data],
        columns=[str(ele) for ele in element_ids],
    )


def _copy_source_components(base_ss):
    """Extract source element IDs, vsource data, and msource data."""
    source_eles = [int(ele) for ele in np.asarray(base_ss.source_eles).reshape(-1)]

    if base_ss.vsource is None:
        if source_eles:
            raise ValueError("base_ss has source elements but no vsource")
        return [], None, []

    vsource_time = _time_array(base_ss.vsource)
    vsource_data = _data_array(base_ss.vsource).copy()

    if vsource_data.shape[1] != len(source_eles):
        raise ValueError(
            f"base_ss source count is {len(source_eles)}, "
            f"but vsource has {vsource_data.shape[1]} columns"
        )

    msource_data_list = []
    for msource in base_ss.msource or []:
        msource_data = _data_array(msource).copy()
        if msource_data.shape[1] != len(source_eles):
            raise ValueError(
                f"base_ss source count is {len(source_eles)}, "
                f"but an msource tracer has {msource_data.shape[1]} columns"
            )
        msource_data_list.append(
            (_time_array(msource), msource_data)
        )

    return source_eles, (vsource_time, vsource_data), msource_data_list


def _copy_sink_components(base_ss):
    """Extract sink element IDs and vsink data."""
    sink_eles_raw = getattr(base_ss, "sink_eles", None)

    if sink_eles_raw is None:
        sink_eles = []
    else:
        sink_eles = [int(ele) for ele in np.asarray(sink_eles_raw).reshape(-1)]

    if base_ss.vsink is None:
        if sink_eles:
            raise ValueError("base_ss has sink elements but no vsink")
        return [], None

    vsink_time = _time_array(base_ss.vsink)
    vsink_data = _data_array(base_ss.vsink).copy()

    if vsink_data.shape[1] != len(sink_eles):
        raise ValueError(
            f"base_ss sink count is {len(sink_eles)}, "
            f"but vsink has {vsink_data.shape[1]} columns"
        )

    return sink_eles, (vsink_time, vsink_data)


def _move_source_element(
    source_eles: list[int],
    old_ele: int,
    target_ele: int,
) -> None:
    """Relabel one source element without changing its time series."""
    if old_ele not in source_eles:
        raise ValueError(f"Source element {old_ele} is not in base_ss")

    if target_ele != old_ele and target_ele in source_eles:
        raise ValueError(
            f"Cannot move source element {old_ele} to {target_ele}: "
            "the target element already contains a source"
        )

    source_eles[source_eles.index(old_ele)] = target_ele


def _move_sink_element(
    sink_eles: list[int],
    old_ele: int,
    target_ele: int,
) -> None:
    """Relabel one sink element without changing its time series."""
    if old_ele not in sink_eles:
        raise ValueError(f"Sink element {old_ele} is not in base_ss")

    if target_ele != old_ele and target_ele in sink_eles:
        raise ValueError(
            f"Cannot move sink element {old_ele} to {target_ele}: "
            "the target element already contains a sink"
        )

    sink_eles[sink_eles.index(old_ele)] = target_ele


def _add_negative_source_part_to_sink(
    target_ele: int,
    source_column_idx: int,
    source_data: np.ndarray,
    source_time: np.ndarray,
    sink_eles: list[int],
    sink_time_and_data: tuple[np.ndarray, np.ndarray] | None,
) -> tuple[list[int], tuple[np.ndarray, np.ndarray] | None, int, float]:
    """
    Move negative source values to a sink at the same element.

    Returns updated sink information, number of negative records, and minimum
    negative flow.
    """
    raw_source = source_data[:, source_column_idx].copy()
    negative_part = np.minimum(raw_source, 0.0)
    n_negative = int(np.sum(negative_part < 0.0))
    min_negative = float(np.min(negative_part)) if len(negative_part) else 0.0

    source_data[:, source_column_idx] = np.maximum(raw_source, 0.0)

    if n_negative == 0:
        return sink_eles, sink_time_and_data, 0, min_negative

    if sink_time_and_data is None:
        sink_time = source_time.copy()
        sink_data = negative_part.reshape(-1, 1)
        sink_eles = [int(target_ele)]
        return sink_eles, (sink_time, sink_data), n_negative, min_negative

    sink_time, sink_data = sink_time_and_data

    if len(sink_time) != len(source_time) or not np.allclose(
        np.asarray(sink_time, dtype=float),
        np.asarray(source_time, dtype=float),
        rtol=0.0,
        atol=1.0e-6,
    ):
        raise ValueError(
            "Cannot add artificial-island negative sink because vsource and "
            "vsink time arrays differ"
        )

    if target_ele in sink_eles:
        sink_idx = sink_eles.index(target_ele)
        sink_data[:, sink_idx] += negative_part
    else:
        sink_eles.append(int(target_ele))
        sink_data = np.column_stack((sink_data, negative_part))

    return sink_eles, (sink_time, sink_data), n_negative, min_negative


def _remove_source_columns(
    source_eles: list[int],
    source_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    msource_data_list: list[tuple[np.ndarray, np.ndarray]],
    remove_eles: set[int],
) -> tuple[
    list[int],
    tuple[np.ndarray, np.ndarray] | None,
    list[tuple[np.ndarray, np.ndarray]],
]:
    """Remove selected source elements from vsource and all msource tracers."""
    if not remove_eles or source_time_and_data is None:
        return source_eles, source_time_and_data, msource_data_list

    keep = np.array(
        [ele not in remove_eles for ele in source_eles],
        dtype=bool,
    )

    source_time, source_data = source_time_and_data
    source_eles = [ele for ele, keep_it in zip(source_eles, keep) if keep_it]
    source_data = source_data[:, keep]

    updated_msource = []
    for tracer_time, tracer_data in msource_data_list:
        updated_msource.append((tracer_time, tracer_data[:, keep]))

    if len(source_eles) == 0:
        return [], None, []

    return (
        source_eles,
        (source_time, source_data),
        updated_msource,
    )


def _remove_sink_columns(
    sink_eles: list[int],
    sink_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    remove_eles: set[int],
) -> tuple[list[int], tuple[np.ndarray, np.ndarray] | None]:
    """Remove selected sink elements from vsink."""
    if not remove_eles or sink_time_and_data is None:
        return sink_eles, sink_time_and_data

    keep = np.array(
        [ele not in remove_eles for ele in sink_eles],
        dtype=bool,
    )

    sink_time, sink_data = sink_time_and_data
    sink_eles = [ele for ele, keep_it in zip(sink_eles, keep) if keep_it]
    sink_data = sink_data[:, keep]

    if len(sink_eles) == 0:
        return [], None

    return sink_eles, (sink_time, sink_data)


def _build_source_sink(
    source_eles: list[int],
    source_time_and_data: tuple[np.ndarray, np.ndarray] | None,
    msource_data_list: list[tuple[np.ndarray, np.ndarray]],
    sink_eles: list[int],
    sink_time_and_data: tuple[np.ndarray, np.ndarray] | None,
) -> source_sink:
    """Build a new source_sink object from patched arrays."""
    if source_time_and_data is None:
        vsource = None
        msource_list = None
    else:
        source_time, source_data = source_time_and_data
        vsource = _build_timehistory(source_time, source_data, source_eles)

        msource_list = []
        for tracer_time, tracer_data in msource_data_list:
            msource_list.append(
                _build_timehistory(tracer_time, tracer_data, source_eles)
            )

    if sink_time_and_data is None:
        vsink = None
    else:
        sink_time, sink_data = sink_time_and_data
        vsink = _build_timehistory(sink_time, sink_data, sink_eles)

    return source_sink(
        vsource=vsource,
        vsink=vsink,
        msource=msource_list,
    )


def _write_diagnostics(
    output_dir: Path,
    patched_ss: source_sink,
    hgrid,
) -> None:
    """Write patched source/sink element coordinates and mean flows."""
    output_dir.mkdir(parents=True, exist_ok=True)
    xctr, yctr = _compute_grid_centers(hgrid)

    source_eles = [
        int(ele) for ele in np.asarray(patched_ss.source_eles).reshape(-1)
    ]
    if patched_ss.vsource is not None and source_eles:
        source_mean = np.mean(_data_array(patched_ss.vsource), axis=0)
        source_xyz = np.c_[
            xctr[np.asarray(source_eles) - 1],
            yctr[np.asarray(source_eles) - 1],
            source_mean,
        ]
        np.savetxt(
            output_dir / "patched_vsource.xyz",
            source_xyz,
            fmt="%.10f %.10f %.8f",
            header="lon lat mean_vsource",
            comments="",
        )

    sink_eles_raw = getattr(patched_ss, "sink_eles", None)
    sink_eles = (
        []
        if sink_eles_raw is None
        else [int(ele) for ele in np.asarray(sink_eles_raw).reshape(-1)]
    )
    if patched_ss.vsink is not None and sink_eles:
        sink_mean = np.mean(_data_array(patched_ss.vsink), axis=0)
        sink_xyz = np.c_[
            xctr[np.asarray(sink_eles) - 1],
            yctr[np.asarray(sink_eles) - 1],
            sink_mean,
        ]
        np.savetxt(
            output_dir / "patched_vsink.xyz",
            sink_xyz,
            fmt="%.10f %.10f %.8f",
            header="lon lat mean_vsink",
            comments="",
        )


def _check_matching_time(
    reference_time: np.ndarray,
    candidate_time: np.ndarray,
    label: str,
) -> None:
    """Require two SCHISM TimeHistory time arrays to match."""
    reference_time = np.asarray(reference_time, dtype=float).reshape(-1)
    candidate_time = np.asarray(candidate_time, dtype=float).reshape(-1)

    if len(reference_time) != len(candidate_time) or not np.allclose(
        reference_time,
        candidate_time,
        rtol=0.0,
        atol=1.0e-6,
    ):
        raise ValueError(
            f"Time array mismatch while restoring {label}: "
            f"base={len(reference_time)} records, "
            f"original={len(candidate_time)} records"
        )



def _as_utc_timestamp(value) -> pd.Timestamp:
    """Convert a datetime-like value to a timezone-aware UTC Timestamp."""
    ts = pd.Timestamp(value)
    if ts.tzinfo is None:
        return ts.tz_localize("UTC")
    return ts.tz_convert("UTC")


def _model_datetimes(start_time, model_time: np.ndarray) -> pd.DatetimeIndex:
    """Convert SCHISM elapsed seconds to UTC datetimes."""
    start = _as_utc_timestamp(start_time)
    seconds = np.asarray(model_time, dtype=float).reshape(-1)
    return pd.DatetimeIndex(start + pd.to_timedelta(seconds, unit="s"))


def _station_id_from_record(record) -> str:
    """Return a normalized USGS station ID from a download record."""
    info = getattr(record, "station_info", {}) or {}
    return str(info.get("id", "")).strip()


def _extract_usgs_values(record) -> pd.Series:
    """
    Return a cleaned, time-indexed USGS value series.

    download_stations() normally returns a dataframe with ``date`` and
    ``value``. A numeric fallback is retained for minor API-format changes.
    """
    df = record.df.copy()
    if "date" not in df.columns:
        raise ValueError("USGS record does not contain a 'date' column")

    index = pd.DatetimeIndex(pd.to_datetime(df["date"], utc=True))

    if "value" in df.columns:
        values = pd.to_numeric(df["value"], errors="coerce")
    else:
        candidate_columns = [
            col for col in df.columns
            if col != "date" and pd.api.types.is_numeric_dtype(
                pd.to_numeric(df[col], errors="coerce")
            )
        ]
        if not candidate_columns:
            raise ValueError("USGS record does not contain a numeric value column")
        values = pd.to_numeric(df[candidate_columns[0]], errors="coerce")

    series = pd.Series(values.to_numpy(dtype=float), index=index)
    series = series.replace([np.inf, -np.inf], np.nan).dropna()
    series = series[~series.index.duplicated(keep="first")].sort_index()
    return series


def _download_usgs_series(
    station_id: str,
    parameter_id: str,
    start_time,
    end_time,
    usgs_cache_folder,
) -> pd.Series | None:
    """Download and validate one USGS station/parameter time series."""
    cache_dir = Path(usgs_cache_folder)
    cache_dir.mkdir(parents=True, exist_ok=True)

    start = _as_utc_timestamp(start_time)
    end = _as_utc_timestamp(end_time)
    padded_start = start - pd.Timedelta(days=1)
    padded_end = end + pd.Timedelta(days=1)

    cache_file = (
        cache_dir
        / (
            f"artificial_island_usgs_{station_id}_{parameter_id}_"
            f"{padded_start.strftime('%Y%m%d')}_"
            f"{padded_end.strftime('%Y%m%d')}.pq"
        )
    )

    try:
        records = download_stations(
            param_id=str(parameter_id),
            station_ids=[str(station_id)],
            cache_fname=str(cache_file),
            datelist=pd.date_range(
                start=padded_start.tz_localize(None),
                end=padded_end.tz_localize(None),
            ),
        )
    except Exception as exc:
        print(
            f"[ARTIFICIAL ISLAND PATCH] warning: USGS download failed "
            f"for station {station_id}, parameter {parameter_id}: {exc}"
        )
        return None

    matching = [
        record
        for record in records
        if _station_id_from_record(record) == str(station_id)
    ]
    if not matching:
        print(
            f"[ARTIFICIAL ISLAND PATCH] warning: USGS station "
            f"{station_id}, parameter {parameter_id} returned no data."
        )
        return None

    record = matching[0]
    try:
        series = _extract_usgs_values(record)
    except Exception as exc:
        print(
            f"[ARTIFICIAL ISLAND PATCH] warning: failed to parse USGS "
            f"station {station_id}, parameter {parameter_id}: {exc}"
        )
        return None

    if series.empty:
        return None

    # Keep partial records. Missing intervals and long gaps are filled later
    # with the original NWM/source forcing rather than rejecting the station.
    return series


def _interpolate_usgs_with_original_fallback(
    series: pd.Series | None,
    target_time: pd.DatetimeIndex,
    original_values: np.ndarray,
    scale: float = 1.0,
    max_interp_gap: pd.Timedelta = MAX_USGS_INTERPOLATION_GAP,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Blend USGS observations with original forcing.

    Exact USGS observations and interpolation across short observation gaps
    are used. Model times outside the observation period, inside long gaps,
    or otherwise unsupported retain original forcing.

    Returns
    -------
    blended
        Final values containing USGS where supported and original forcing
        elsewhere.
    use_usgs
        Boolean mask identifying records supplied by USGS.
    """
    original_values = np.asarray(original_values, dtype=float).reshape(-1)
    target_time = pd.DatetimeIndex(target_time)

    if len(original_values) != len(target_time):
        raise ValueError(
            "Original forcing length does not match model datetime length"
        )

    blended = original_values.copy()
    use_usgs = np.zeros(len(target_time), dtype=bool)

    if series is None:
        return blended, use_usgs

    series = series.copy()
    series = series.replace([np.inf, -np.inf], np.nan).dropna()
    series = series[~series.index.duplicated(keep="first")].sort_index()

    if series.empty:
        return blended, use_usgs

    if series.index.tz is None:
        series.index = series.index.tz_localize("UTC")
    else:
        series.index = series.index.tz_convert("UTC")

    if target_time.tz is None:
        target_time = target_time.tz_localize("UTC")
    else:
        target_time = target_time.tz_convert("UTC")

    obs_time_ns = series.index.asi8
    target_time_ns = target_time.asi8
    obs_values = series.to_numpy(dtype=float) * float(scale)

    right_idx = np.searchsorted(
        obs_time_ns,
        target_time_ns,
        side="left",
    )
    left_idx = right_idx - 1

    safe_right = np.minimum(right_idx, len(obs_time_ns) - 1)
    exact = (
        (right_idx < len(obs_time_ns))
        & (obs_time_ns[safe_right] == target_time_ns)
    )

    if np.any(exact):
        exact_obs_idx = right_idx[exact]
        blended[exact] = obs_values[exact_obs_idx]
        use_usgs[exact] = True

    between = (
        (left_idx >= 0)
        & (right_idx < len(obs_time_ns))
        & (~exact)
    )

    candidate_idx = np.where(between)[0]
    if len(candidate_idx) > 0:
        left_obs_idx = left_idx[candidate_idx]
        right_obs_idx = right_idx[candidate_idx]

        gap_ns = (
            obs_time_ns[right_obs_idx]
            - obs_time_ns[left_obs_idx]
        )
        short_gap = gap_ns <= int(max_interp_gap.value)
        valid_target_idx = candidate_idx[short_gap]

        if len(valid_target_idx) > 0:
            left_obs_idx = left_idx[valid_target_idx]
            right_obs_idx = right_idx[valid_target_idx]

            denominator = (
                obs_time_ns[right_obs_idx]
                - obs_time_ns[left_obs_idx]
            )
            fraction = (
                target_time_ns[valid_target_idx]
                - obs_time_ns[left_obs_idx]
            ) / denominator

            interpolated = (
                obs_values[left_obs_idx]
                + fraction
                * (
                    obs_values[right_obs_idx]
                    - obs_values[left_obs_idx]
                )
            )

            blended[valid_target_idx] = interpolated
            use_usgs[valid_target_idx] = True

    return blended, use_usgs


def _get_usgs_forcing(
    name: str,
    start_time,
    model_time: np.ndarray,
    usgs_cache_folder,
) -> tuple[pd.Series | None, pd.Series | None, str | None]:
    """
    Return raw USGS flow and temperature series plus station ID.

    The caller blends these series with original forcing on the model time
    axis, retaining original NWM/source values outside valid USGS coverage.
    """
    station_id = USGS_STATION_BY_NAME.get(name)
    if not station_id:
        print(
            f"[ARTIFICIAL ISLAND PATCH] warning: {name} has "
            "use_usgs_obs=true but no station ID in USGS_STATION_BY_NAME; "
            "original forcing will be used."
        )
        return None, None, None

    target_time = _model_datetimes(start_time, model_time)
    period_start = target_time[0]
    period_end = target_time[-1]

    flow_series = _download_usgs_series(
        station_id=station_id,
        parameter_id=USGS_FLOW_PARAMETER_ID,
        start_time=period_start,
        end_time=period_end,
        usgs_cache_folder=usgs_cache_folder,
    )

    temperature_series = _download_usgs_series(
        station_id=station_id,
        parameter_id=USGS_TEMPERATURE_PARAMETER_ID,
        start_time=period_start,
        end_time=period_end,
        usgs_cache_folder=usgs_cache_folder,
    )

    return flow_series, temperature_series, station_id




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
    manual_nwm2usgs = {
        19406836: "07381490",
        15708755: "02489500",
        18928090: "07375175",
        19269176: "07374000",
        16665157: "02244040",
        2590217: "01463500",
        6186156: "01358000",
    }
    for feature_id, station_id in manual_nwm2usgs.items():
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
    """Download USGS 00010 for all required stations and return usable series."""
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

    result: dict[str, pd.Series] = {}
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
) -> tuple[
    tuple[np.ndarray, np.ndarray] | None,
    list[tuple[np.ndarray, np.ndarray]],
    list[int],
    tuple[np.ndarray, np.ndarray] | None,
]:
    """
    Replace flow and/or temperature in one existing relocated source column.

    The source element and source-column count are unchanged. USGS records are
    blended with the existing relocated forcing: unsupported times retain the
    existing values.
    """
    name = point["name"]
    radius_m = point["max_search_radius_m"]

    if source_time_and_data is None or not source_eles:
        raise ValueError(
            f"[ARTIFICIAL ISLAND PATCH] {name}: relocated base has no source"
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
            f"[ARTIFICIAL ISLAND PATCH] {name}: no relocated source found "
            f"within {radius_m:.1f} m; nearest={nearest_text}. "
            "Replace-only entries never create a new source."
        )

    station_id = USGS_STATION_BY_NAME.get(name)
    if not station_id:
        raise ValueError(
            f"[ARTIFICIAL ISLAND PATCH] {name}: no USGS station is configured"
        )

    source_time, source_data = source_time_and_data
    source_idx = source_eles.index(target_ele)
    target_datetime = _model_datetimes(start_time, source_time)

    flow_series, temperature_series, _ = _get_usgs_forcing(
        name=name,
        start_time=start_time,
        model_time=source_time,
        usgs_cache_folder=usgs_cache_folder,
    )

    if point["replace_flow"]:
        if flow_series is None:
            print(
                f"[ARTIFICIAL ISLAND PATCH] warning: {name}: USGS station "
                f"{station_id} returned no flow; relocated vsource is unchanged."
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
                f"[ARTIFICIAL ISLAND PATCH] {name}: replaced relocated "
                f"vsource at element {target_ele} using station {station_id} "
                f"for {int(use_usgs_flow.sum())}/{len(use_usgs_flow)} records; "
                f"existing relocated forcing retained for "
                f"{int((~use_usgs_flow).sum())} records."
            )

    if point["replace_temperature"]:
        if not msource_data_list:
            print(
                f"[ARTIFICIAL ISLAND PATCH] warning: {name}: base_ss has no "
                "msource tracer; temperature replacement was skipped."
            )
        elif temperature_series is None:
            print(
                f"[ARTIFICIAL ISLAND PATCH] warning: {name}: USGS station "
                f"{station_id} returned no temperature; relocated temperature "
                "msource is unchanged."
            )
        else:
            tracer_time, tracer_data = msource_data_list[0]
            _check_matching_time(
                source_time,
                tracer_time,
                f"relocated temperature replacement {name}",
            )

            original_temperature = tracer_data[:, source_idx].copy()
            replaced_temperature, use_usgs_temperature = (
                _interpolate_usgs_with_original_fallback(
                    series=temperature_series,
                    target_time=target_datetime,
                    original_values=original_temperature,
                    scale=1.0,
                )
            )
            tracer_data[:, source_idx] = replaced_temperature
            msource_data_list[0] = (tracer_time, tracer_data)

            print(
                f"[ARTIFICIAL ISLAND PATCH] {name}: replaced relocated "
                f"temperature msource tracer 1 at element {target_ele} using "
                f"station {station_id} for "
                f"{int(use_usgs_temperature.sum())}/"
                f"{len(use_usgs_temperature)} records; existing temperature "
                f"retained for {int((~use_usgs_temperature).sum())} records."
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
                f"[ARTIFICIAL ISLAND PATCH] {name}: moved "
                f"{n_negative} negative relocated-source record(s) to vsink "
                f"at element {target_ele}; minimum={min_negative:.6f} m3/s."
            )
        else:
            print(
                f"[ARTIFICIAL ISLAND PATCH] {name}: "
                "allow_negative_sink enabled; no negative records."
            )

    print(
        f"[ARTIFICIAL ISLAND PATCH] {name}: replace-only operation matched "
        f"relocated source element {target_ele} at {distance_m:.1f} m."
    )

    return (
        (source_time, source_data),
        msource_data_list,
        sink_eles,
        sink_time_and_data,
    )

def _append_source_column(
    source_eles: list[int],
    source_time_and_data,
    msource_data_list,
    target_ele: int,
    source_time: np.ndarray,
    flow: np.ndarray,
    tracer_columns: list[np.ndarray],
):
    """Append one new source and corresponding msource tracer columns."""
    if target_ele in source_eles:
        raise ValueError(
            f"Target element {target_ele} already contains a source"
        )

    flow = np.asarray(flow, dtype=float).reshape(-1)
    source_time = np.asarray(source_time, dtype=float).reshape(-1)

    if len(flow) != len(source_time):
        raise ValueError("New source flow length does not match model time")

    if source_time_and_data is None:
        combined_source_time = source_time.copy()
        combined_source_data = flow.reshape(-1, 1)
    else:
        combined_source_time, combined_source_data = source_time_and_data
        _check_matching_time(
            combined_source_time,
            source_time,
            f"new source element {target_ele}",
        )
        combined_source_data = np.column_stack(
            (combined_source_data, flow)
        )

    if len(msource_data_list) != len(tracer_columns):
        raise ValueError(
            "New source tracer count does not match base msource tracer count"
        )

    updated_msource = []
    for tracer_idx, (
        (tracer_time, tracer_data),
        tracer_column,
    ) in enumerate(zip(msource_data_list, tracer_columns)):
        _check_matching_time(
            tracer_time,
            source_time,
            f"new source msource tracer {tracer_idx + 1}",
        )
        tracer_column = np.asarray(tracer_column, dtype=float).reshape(-1)
        if len(tracer_column) != len(source_time):
            raise ValueError(
                f"Tracer {tracer_idx + 1} length does not match model time"
            )
        updated_msource.append(
            (
                tracer_time,
                np.column_stack((tracer_data, tracer_column)),
            )
        )

    source_eles.append(int(target_ele))
    return (
        source_eles,
        (combined_source_time, combined_source_data),
        updated_msource,
    )


def _append_sink_column(
    sink_eles: list[int],
    sink_time_and_data,
    target_ele: int,
    sink_time: np.ndarray,
    sink_values: np.ndarray,
):
    """Append one new sink column."""
    if target_ele in sink_eles:
        raise ValueError(
            f"Target element {target_ele} already contains a sink"
        )

    sink_time = np.asarray(sink_time, dtype=float).reshape(-1)
    sink_values = np.asarray(sink_values, dtype=float).reshape(-1)

    if len(sink_time) != len(sink_values):
        raise ValueError("New sink length does not match model time")

    if sink_time_and_data is None:
        combined_time = sink_time.copy()
        combined_data = sink_values.reshape(-1, 1)
    else:
        combined_time, combined_data = sink_time_and_data
        _check_matching_time(
            combined_time,
            sink_time,
            f"new sink element {target_ele}",
        )
        combined_data = np.column_stack(
            (combined_data, sink_values)
        )

    sink_eles.append(int(target_ele))
    return sink_eles, (combined_time, combined_data)



def patch_artificial_island_source_sink(
    base_ss: source_sink,
    hgrid,
    original_source_sink_dir,
    patch_info_file: str | Path | dict,
    start_time=None,
    rnday=None,
    usgs_cache_folder=None,
    output_dir: str | Path | None = None,
    relocated_source_sink_dir: str | Path | None = None,
    nwm_shapefile: str | Path | None = None,
    states=None,
    replace_all_source_temperature: bool = True,
) -> source_sink:
    """
    Restore only YAML-listed artificial-island source/sink forcing.

    The forcing columns are copied from ``original_source_sink_dir`` and
    appended to the already relocated ``base_ss``. Existing relocated source
    or sink columns are never moved.

    Before YAML-specific operations, all relocated sources can receive
    automatic USGS temperature replacement using relocated sources.json and
    the same upstream NWM/USGS search logic as source_nwm2usgs(). Explicit
    USGS_STATION_BY_NAME replacements are applied afterward as overrides.
    ``rnday`` is retained for workflow compatibility.
    """
    if base_ss is None:
        raise ValueError("base_ss must not be None")

    _ = rnday
    if start_time is None:
        raise ValueError("start_time is required for artificial-island patching")
    usgs_cache_folder = Path(
        usgs_cache_folder or Path(output_dir or Path.cwd()) / "USGS_cache"
    )

    original_source_sink_dir = Path(original_source_sink_dir)
    relocated_source_sink_dir = Path(
        relocated_source_sink_dir
        or original_source_sink_dir.parent / "relocated_source_sink"
    )
    nwm_shapefile = Path(nwm_shapefile)
    states = list(states or STOFS3D_ATL_STATES)

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
    force_points = _normalize_force_points(patch_info)
    replace_relocated_points = _normalize_replace_relocated_points(
        patch_info
    )
    large_constant_sink_points = (
        _normalize_large_constant_sink_points(patch_info)
    )
    exclude_points = _normalize_exclude_points(patch_info)

    if (
        not force_points
        and not replace_relocated_points
        and not large_constant_sink_points
        and not exclude_points
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
    replaced_relocated_count = 0
    large_constant_sink_count = 0
    automatically_temperature_replaced_count = 0

    # ------------------------------------------------------------------
    # First, apply automatic USGS temperature replacement to every
    # relocated source using the same upstream NWM/USGS association logic
    # as source_nwm2usgs(). Explicit YAML replacements are applied later
    # and therefore override this automatic result.
    # ------------------------------------------------------------------
    if replace_all_source_temperature:
        diagnostics_dir = (
            None
            if output_dir is None
            else Path(output_dir) / "automatic_temperature"
        )
        if diagnostics_dir is not None:
            diagnostics_dir.mkdir(parents=True, exist_ok=True)

        (
            msource_data_list,
            automatically_temperature_replaced_count,
        ) = _replace_all_relocated_source_temperatures(
            source_eles=source_eles,
            source_time_and_data=source_time_and_data,
            msource_data_list=msource_data_list,
            xctr=xctr,
            yctr=yctr,
            relocated_source_sink_dir=relocated_source_sink_dir,
            start_time=start_time,
            usgs_cache_folder=usgs_cache_folder,
            nwm_shapefile=nwm_shapefile,
            states=states,
            diagnostics_dir=diagnostics_dir,
        )

    # ------------------------------------------------------------------
    # Replace values in existing relocated sources.
    # No source element is added, removed, or moved in this step.
    # ------------------------------------------------------------------
    for point in replace_relocated_points:
        (
            source_time_and_data,
            msource_data_list,
            sink_eles,
            sink_time_and_data,
        ) = _replace_existing_relocated_source(
            point=point,
            source_eles=source_eles,
            source_time_and_data=source_time_and_data,
            msource_data_list=msource_data_list,
            sink_eles=sink_eles,
            sink_time_and_data=sink_time_and_data,
            xctr=xctr,
            yctr=yctr,
            transformer=transformer,
            start_time=start_time,
            usgs_cache_folder=usgs_cache_folder,
        )
        replaced_relocated_count += 1

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
            station_id = None

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
                    station_id,
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
                    f"[ARTIFICIAL ISLAND PATCH] {name}: USGS station "
                    f"{station_id} supplied vsource for "
                    f"{int(use_usgs_flow.sum())}/{len(use_usgs_flow)} "
                    "records; original NWM supplied "
                    f"{int((~use_usgs_flow).sum())} records."
                )

                if tracer_columns:
                    print(
                        f"[ARTIFICIAL ISLAND PATCH] {name}: USGS station "
                        f"{station_id} supplied temperature msource for "
                        f"{int(use_usgs_temperature.sum())}/"
                        f"{len(use_usgs_temperature)} records; original "
                        "msource supplied "
                        f"{int((~use_usgs_temperature).sum())} records."
                    )

                source_origin_message = (
                    f"USGS station {station_id} blended with original "
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

    # ------------------------------------------------------------------
    # Remove source/sink forcing near YAML exclusion points.
    # ------------------------------------------------------------------
    remove_sources: set[int] = set()
    remove_sinks: set[int] = set()

    for point in exclude_points:
        if "source" in point["remove"]:
            found = _elements_within_radius(
                x=point["x"],
                y=point["y"],
                radius_m=point["radius_m"],
                candidate_element_ids=source_eles,
                xctr=xctr,
                yctr=yctr,
                transformer=transformer,
            )
            remove_sources.update(found)
            print(
                f"[ARTIFICIAL ISLAND PATCH] {point['name']}: "
                f"marked {len(found)} source(s) for removal: {found}"
            )

        if "sink" in point["remove"]:
            found = _elements_within_radius(
                x=point["x"],
                y=point["y"],
                radius_m=point["radius_m"],
                candidate_element_ids=sink_eles,
                xctr=xctr,
                yctr=yctr,
                transformer=transformer,
            )
            remove_sinks.update(found)
            print(
                f"[ARTIFICIAL ISLAND PATCH] {point['name']}: "
                f"marked {len(found)} sink(s) for removal: {found}"
            )

    (
        source_eles,
        source_time_and_data,
        msource_data_list,
    ) = _remove_source_columns(
        source_eles=source_eles,
        source_time_and_data=source_time_and_data,
        msource_data_list=msource_data_list,
        remove_eles=remove_sources,
    )

    sink_eles, sink_time_and_data = _remove_sink_columns(
        sink_eles=sink_eles,
        sink_time_and_data=sink_time_and_data,
        remove_eles=remove_sinks,
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
        "[ARTIFICIAL ISLAND PATCH] automatically temperature-replaced "
        f"{automatically_temperature_replaced_count} relocated source(s)."
    )
    print(
        "[ARTIFICIAL ISLAND PATCH] replaced "
        f"{replaced_relocated_count} existing relocated source(s)."
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

def remove_overlapping_background_sinks(
    background_ss: source_sink,
    base_ss: source_sink,
    output_dir: str | Path | None = None,
) -> source_sink:
    """
    Remove background-sink columns whose elements are already used by base_ss.

    This preserves the original set_constant_sink() function and applies the
    overlap removal afterward.

    Both base source elements and base sink elements are excluded from the
    background sink object.
    """
    if background_ss is None:
        raise ValueError("background_ss must not be None")
    if base_ss is None:
        raise ValueError("base_ss must not be None")

    base_source_eles = {
        int(ele) for ele in np.asarray(base_ss.source_eles).reshape(-1)
    }

    base_sink_raw = getattr(base_ss, "sink_eles", None)
    base_sink_eles = (
        set()
        if base_sink_raw is None
        else {int(ele) for ele in np.asarray(base_sink_raw).reshape(-1)}
    )

    excluded = base_source_eles | base_sink_eles

    background_sink_raw = getattr(background_ss, "sink_eles", None)
    background_sink_eles = (
        []
        if background_sink_raw is None
        else [
            int(ele)
            for ele in np.asarray(background_sink_raw).reshape(-1)
        ]
    )

    if background_ss.vsink is None or not background_sink_eles:
        return background_ss

    sink_time = _time_array(background_ss.vsink)
    sink_data = _data_array(background_ss.vsink)

    keep = np.array(
        [ele not in excluded for ele in background_sink_eles],
        dtype=bool,
    )

    kept_eles = [
        ele
        for ele, keep_it in zip(background_sink_eles, keep)
        if keep_it
    ]
    kept_data = sink_data[:, keep]

    removed = [
        ele
        for ele, keep_it in zip(background_sink_eles, keep)
        if not keep_it
    ]

    print(
        "[BACKGROUND SINK PATCH] removed "
        f"{len(removed)} overlapping background sink element(s)."
    )

    patched_background_ss = source_sink(
        vsource=None,
        vsink=_build_timehistory(sink_time, kept_data, kept_eles),
        msource=None,
    )

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        patched_background_ss.writer(str(output_dir))

    return patched_background_ss


__all__ = [
    "patch_artificial_island_source_sink",
    "zero_artificial_island_sources_after_replace_USGS_before_relocation",
    "remove_overlapping_background_sinks",
]
