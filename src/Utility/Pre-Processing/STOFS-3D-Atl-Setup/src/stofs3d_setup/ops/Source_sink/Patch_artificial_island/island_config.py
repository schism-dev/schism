"""Normalize configuration owned by artificial-island corrections."""

from stofs3d_setup.ops.Source_sink.correction_config import (
    _as_dict,
    load_source_sink_corrections,
)


def _load_patch_info(patch_info_file):
    """Load artificial-island settings from YAML or a mapping."""
    return load_source_sink_corrections(patch_info_file)


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

        use_usgs_obs = bool(p.get("use_usgs_obs", False))
        flow_station_id = p.get("flow_station_id")
        temperature_station_id = p.get("temperature_station_id")
        for key, station_id in (
            ("flow_station_id", flow_station_id),
            ("temperature_station_id", temperature_station_id),
        ):
            if station_id is not None and (
                not isinstance(station_id, str) or not station_id.strip()
            ):
                raise ValueError(
                    f"{key} for artificial-island entry "
                    f"{p.get('name', '')!r} must be a quoted station ID "
                    "or null"
                )
        if use_usgs_obs and not flow_station_id:
            raise ValueError(
                f"Artificial-island entry {p.get('name', '')!r} enables "
                "USGS observations but has no flow_station_id"
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
                "use_usgs_obs": use_usgs_obs,
                "flow_station_id": flow_station_id,
                "temperature_station_id": temperature_station_id,
            }
        )

    return points


def _normalize_large_constant_sink_points(patch_info: dict) -> list[dict]:
    """Normalize artificial-island locations for large constant sinks."""
    raw_points = patch_info.get(
        "large_constant_sink_artificial_island_locations"
    ) or []

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
