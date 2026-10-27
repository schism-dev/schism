"""Direct USGS flow and temperature hard resets for existing sources."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer
from shapely.geometry import Point
from shapely.ops import transform as transform_geometry

from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    CFS_TO_CMS,
    NWM_TO_USGS_FLOW_ADJUSTMENT,
    NWM_TO_USGS_TEMPERATURE_SEARCH,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.nwm_usgs_station_search import (
    associate_poi_with_nwm,
    preprocess_nwm_shp,
)
from stofs3d_setup.ops.USGS.ogc.download_usgs import (
    discover_monitoring_locations,
    normalize_station_id,
)
from stofs3d_setup.utils.utils import STOFS3D_ATL_STATES
from stofs3d_setup.ops.Source_sink.workflow_constants import (
    NWM_ECGC_SHAPEFILE,
)
from stofs3d_setup.ops.USGS.ogc.series import (
    _get_usgs_flow,
    _get_usgs_temperature,
    _interpolate_usgs_with_original_fallback,
    _model_datetimes,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.source_temperature import (
    FIRST_USABLE,
    _blend_complete_temperature_or_preserve,
    _mixed_temperature_columns,
    _validate_temperature_pooling,
)
from stofs3d_setup.ops.Source_sink.correction_config import load_source_sink_corrections
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.hard_reset_config import (
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


def _load_source_feature_mapping(source_mapping_dir) -> dict[int, list[int]]:
    """Load current source-element to NWM-feature mappings."""
    mapping_file = Path(source_mapping_dir) / "sources.json"
    if not mapping_file.is_file():
        raise FileNotFoundError(
            f"sources.json is required for downstream USGS hard resets: "
            f"{mapping_file}"
        )

    with mapping_file.open("r", encoding="utf-8") as source_file:
        raw_mapping = json.load(source_file)

    mapping = {}
    for element_id, feature_ids in raw_mapping.items():
        if not isinstance(feature_ids, list):
            feature_ids = [feature_ids]
        mapping[int(element_id)] = [
            int(feature_id) for feature_id in feature_ids
        ]
    return mapping


def _station_feature_id(point: dict, nwm_network) -> tuple[str, int]:
    """Associate the configured USGS station with its NWM segment."""
    if point.get("replace_flow", True):
        station_id = point.get("flow_station_id")
        station_links = NWM_TO_USGS_FLOW_ADJUSTMENT
    else:
        station_id = point.get("temperature_station_id")
        station_links = NWM_TO_USGS_TEMPERATURE_SEARCH

    if not station_id:
        station_id = point.get("temperature_station_id")
        station_links = NWM_TO_USGS_TEMPERATURE_SEARCH
    if not station_id:
        raise ValueError(
            f"[SOURCE OVERRIDE] {point['name']}: no station is configured "
            "to locate its NWM route"
        )

    station_id = normalize_station_id(station_id)
    manual_features = [
        int(feature_id)
        for feature_id, linked_station_id in station_links.items()
        if normalize_station_id(linked_station_id) == station_id
    ]
    if len(manual_features) > 1:
        raise ValueError(
            f"[SOURCE OVERRIDE] station {station_id} maps to multiple NWM "
            f"features: {manual_features}"
        )
    if manual_features:
        feature_id = manual_features[0]
        if feature_id not in nwm_network.index:
            raise ValueError(
                f"[SOURCE OVERRIDE] station {station_id} maps to NWM feature "
                f"{feature_id}, which is absent from the NWM network"
            )
        return station_id, feature_id

    locations = discover_monitoring_locations(states=STOFS3D_ATL_STATES)
    station_ids = locations["station_id"].map(normalize_station_id)
    matches = locations.loc[station_ids == station_id]
    if matches.empty:
        raise ValueError(
            f"[SOURCE OVERRIDE] USGS station {station_id} has no location "
            "in the Atlantic station inventory"
        )

    station = matches.iloc[0]
    station_network = associate_poi_with_nwm(
        nwm_network,
        poi=np.asarray(
            [[float(station["lon"]), float(station["lat"])]],
        ),
        poi_names=[station_id],
        poi_label="hard_reset_station",
        invalid_poi_output=None,
    )
    matched_features = station_network.index[
        station_network["hard_reset_station"] == station_id
    ].to_list()
    if len(matched_features) != 1:
        raise ValueError(
            f"[SOURCE OVERRIDE] USGS station {station_id} could not be "
            "associated with exactly one NWM segment"
        )
    return station_id, int(matched_features[0])


def _search_downstream_source(
    start_feature_id: int,
    nwm_network,
    source_feature_mapping: dict[int, list[int]],
    source_element_ids: set[int],
    xctr: np.ndarray,
    yctr: np.ndarray,
    transformer: Transformer,
    max_search_radius_m: float,
) -> tuple[int, int, int, float]:
    """Follow NWM downstream to the first mapped in-grid source element."""
    network_transformer = Transformer.from_crs(
        nwm_network.crs or "EPSG:4326",
        "EPSG:3857",
        always_xy=True,
    )
    visited = set()
    current_feature_id = int(start_feature_id)
    while current_feature_id in nwm_network.index:
        if current_feature_id in visited:
            raise ValueError(
                "[SOURCE OVERRIDE] NWM downstream route contains a cycle "
                f"at feature {current_feature_id}"
            )
        visited.add(current_feature_id)

        mapped_elements = [
            element_id
            for element_id in source_feature_mapping.get(
                current_feature_id, []
            )
            if element_id in source_element_ids
            and 0 < element_id <= len(xctr)
        ]
        if mapped_elements:
            flowline = transform_geometry(
                network_transformer.transform,
                nwm_network.loc[current_feature_id, "geometry"],
            )
            source_x, source_y = transformer.transform(
                xctr[np.asarray(mapped_elements, dtype=int) - 1],
                yctr[np.asarray(mapped_elements, dtype=int) - 1],
            )
            distances = [
                float(flowline.distance(Point(float(x), float(y))))
                for x, y in zip(source_x, source_y)
            ]
            nearest_index = int(np.argmin(distances))
            distance_m = distances[nearest_index]
            if distance_m > max_search_radius_m:
                raise ValueError(
                    "[SOURCE OVERRIDE] first downstream grid source on "
                    f"feature {current_feature_id} is "
                    f"{distance_m:.1f} m from its NWM flowline, outside the "
                    f"{max_search_radius_m:.1f} m search radius"
                )
            return (
                mapped_elements[nearest_index],
                current_feature_id,
                len(visited),
                distance_m,
            )

        downstream = nwm_network.loc[current_feature_id, "to"]
        if pd.isna(downstream):
            break
        next_feature_id = int(downstream)
        if next_feature_id == current_feature_id:
            break
        current_feature_id = next_feature_id

    raise ValueError(
        "[SOURCE OVERRIDE] no downstream source was found from NWM feature "
        f"{start_feature_id} after visiting {len(visited)} segment(s)"
    )


def resolve_hard_reset_source_targets(
    base_ss,
    hgrid,
    points: list[dict],
    source_mapping_dir,
    nwm_shapefile=NWM_ECGC_SHAPEFILE,
) -> list[dict]:
    """Resolve hard resets to the first matching source downstream of USGS."""
    if not points:
        return points

    nwm_network = preprocess_nwm_shp(str(nwm_shapefile))
    source_mapping = _load_source_feature_mapping(source_mapping_dir)
    source_element_ids = {
        int(element_id) for element_id in base_ss.source_eles
    }
    feature_to_elements: dict[int, list[int]] = {}
    for element_id, feature_ids in source_mapping.items():
        for feature_id in feature_ids:
            feature_to_elements.setdefault(feature_id, []).append(element_id)

    xctr, yctr = _compute_grid_centers(hgrid)
    transformer = _make_transformer()
    resolved_points = []
    for point in points:
        station_id, station_feature_id = _station_feature_id(
            point, nwm_network
        )
        (
            target_element_id,
            boundary_feature_id,
            route_length,
            distance_m,
        ) = _search_downstream_source(
            start_feature_id=station_feature_id,
            nwm_network=nwm_network,
            source_feature_mapping=feature_to_elements,
            source_element_ids=source_element_ids,
            xctr=xctr,
            yctr=yctr,
            transformer=transformer,
            max_search_radius_m=point["max_search_radius_m"],
        )
        print(
            f"[SOURCE OVERRIDE] {point['name']}: station {station_id} "
            f"is on NWM feature {station_feature_id}; downstream grid "
            f"source is element {target_element_id} on feature "
            f"{boundary_feature_id} after {route_length} segment(s), "
            f"{distance_m:.1f} m from its flowline."
        )
        resolved_points.append(
            {
                **point,
                "_station_id": station_id,
                "_station_nwm_feature_id": station_feature_id,
                "_target_element_id": target_element_id,
                "_target_distance_m": distance_m,
                "_target_nwm_feature_id": boundary_feature_id,
                "_route_segment_count": route_length,
            }
        )

    route_manifest = [
        {
            "name": point["name"],
            "station_id": point["_station_id"],
            "station_nwm_feature_id": point[
                "_station_nwm_feature_id"
            ],
            "search_direction": "downstream",
            "target_nwm_feature_id": point["_target_nwm_feature_id"],
            "target_source_element_id": point["_target_element_id"],
            "route_segment_count": point["_route_segment_count"],
            "source_to_flowline_distance_m": point[
                "_target_distance_m"
            ],
            "max_search_radius_m": point["max_search_radius_m"],
        }
        for point in resolved_points
    ]
    (Path(source_mapping_dir) / "hard_reset_source_routes.json").write_text(
        json.dumps(route_manifest, indent=2),
        encoding="utf-8",
    )
    return resolved_points


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
    A temperature override is atomic in either pooling mode: it is applied
    only if every model time receives a finite numeric temperature.
    """
    name = point["name"]
    radius_m = point["max_search_radius_m"]
    temperature_pooling = _validate_temperature_pooling(temperature_pooling)

    if source_time_and_data is None or not source_eles:
        raise ValueError(
            f"[SOURCE OVERRIDE] {name}: relocated base has no source"
        )

    target_ele = point.get("_target_element_id")
    if target_ele is None:
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
    else:
        target_ele = int(target_ele)
        if target_ele not in source_eles:
            raise ValueError(
                f"[SOURCE OVERRIDE] {name}: downstream NWM target element "
                f"{target_ele} is not an existing relocated source"
            )
        distance_m = float(point["_target_distance_m"])
        if distance_m > radius_m:
            raise ValueError(
                f"[SOURCE OVERRIDE] {name}: downstream source element "
                f"{target_ele} is {distance_m:.1f} m from the NWM flowline, "
                f"outside the {radius_m:.1f} m search radius"
            )

    flow_station_id = point.get("flow_station_id")
    temperature_station_id = point.get("temperature_station_id")

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
    if point["replace_flow"]:
        flow_series, flow_station_id = _get_usgs_flow(
            name=name,
            station_id=flow_station_id,
            start_time=start_time,
            model_time=source_time,
            usgs_cache_folder=usgs_cache_folder,
        )
    if point["replace_temperature"]:
        temperature_series, temperature_station_id = _get_usgs_temperature(
            name=name,
            station_id=temperature_station_id,
            start_time=start_time,
            model_time=source_time,
            usgs_cache_folder=usgs_cache_folder,
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
            (
                replaced_temperature,
                use_usgs_temperature,
                temperature_applied,
            ) = _blend_complete_temperature_or_preserve(
                series=temperature_series,
                target_time=target_datetime,
                original_values=original_temperature,
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

    target_description = (
        f"NWM feature {point['_target_nwm_feature_id']}"
        if "_target_nwm_feature_id" in point
        else "configured point"
    )
    print(
        f"[SOURCE OVERRIDE] {name}: replace-only operation matched "
        f"relocated source element {target_ele} via {target_description}; "
        f"distance={distance_m:.1f} m."
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
    targets_resolved = all(
        point.get("_target_element_id") is not None for point in points
    )
    if targets_resolved:
        xctr = np.empty(0)
        yctr = np.empty(0)
    else:
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

    if msource_values:
        mixed_columns = np.flatnonzero(
            _mixed_temperature_columns(msource_values[0][1])
        )
        if len(mixed_columns) > 0:
            mixed_elements = [source_eles[idx] for idx in mixed_columns]
            raise ValueError(
                "Source-temperature hard resets left columns that are not "
                "uniformly all-ambient or "
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


def apply_hard_reset_sources_with_usgs(
    base_ss,
    hgrid,
    enabled: bool,
    correction_info,
    start_time,
    usgs_cache_folder,
    temperature_pooling: str = FIRST_USABLE,
    source_mapping_dir=None,
    nwm_shapefile=NWM_ECGC_SHAPEFILE,
):
    """Apply configured direct USGS flow and temperature hard resets.

    This wrapper is available to workflows that explicitly opt into the
    hard-reset corrections. When source mappings are available, each station
    is matched to its first in-grid source along the downstream NWM route.
    """
    if not enabled:
        return base_ss, 0, 0
    if correction_info is None:
        raise ValueError(
            "correction_info is required when USGS hard resets "
            "are enabled"
        )

    points = load_source_override_points(correction_info)
    if source_mapping_dir is not None:
        points = resolve_hard_reset_source_targets(
            base_ss=base_ss,
            hgrid=hgrid,
            points=points,
            source_mapping_dir=source_mapping_dir,
            nwm_shapefile=nwm_shapefile,
        )
    base_ss, flow_count = apply_source_flow_overrides(
        base_ss=base_ss,
        hgrid=hgrid,
        points=points,
        start_time=start_time,
        usgs_cache_folder=usgs_cache_folder,
    )
    base_ss, temperature_count = apply_source_temperature_overrides(
        base_ss=base_ss,
        hgrid=hgrid,
        points=points,
        start_time=start_time,
        usgs_cache_folder=usgs_cache_folder,
        temperature_pooling=temperature_pooling,
    )
    return base_ss, flow_count, temperature_count


def load_source_override_points(correction_info) -> list[dict]:
    """Load and normalize configured source hard-reset locations."""
    corrections = load_source_sink_corrections(correction_info)
    return source_override_points(corrections)
