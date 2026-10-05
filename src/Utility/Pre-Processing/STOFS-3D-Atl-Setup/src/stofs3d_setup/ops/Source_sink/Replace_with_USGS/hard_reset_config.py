"""Normalize configured USGS hard resets for existing sources."""

from stofs3d_setup.ops.Source_sink.correction_config import _as_dict


def _normalize_replace_relocated_points(patch_info: dict) -> list[dict]:
    """Normalize configured existing-source hard-reset locations."""
    raw_points = (
        patch_info.get("hard_reset_source_locations")
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
        replace_flow = bool(p.get("replace_flow", True))
        replace_temperature = bool(p.get("replace_temperature", True))
        flow_station_id = p.get("flow_station_id")
        temperature_station_id = p.get("temperature_station_id")

        if replace_flow and not _valid_station_id(flow_station_id):
            raise ValueError(
                f"{name!r} enables flow hard reset but has no valid "
                "flow_station_id"
            )
        if replace_temperature and not _valid_station_id(
            temperature_station_id
        ):
            raise ValueError(
                f"{name!r} enables temperature hard reset but has no valid "
                "temperature_station_id"
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
                "flow_station_id": flow_station_id,
                "temperature_station_id": temperature_station_id,
                "replace_flow": replace_flow,
                "replace_temperature": replace_temperature,
                "allow_negative_sink": bool(
                    p.get("allow_negative_sink", False)
                ),
            }
        )

    return points


def _valid_station_id(station_id) -> bool:
    """Return whether a station ID is a non-empty string."""
    return isinstance(station_id, str) and bool(station_id.strip())


def source_override_points(corrections: dict) -> list[dict]:
    """Return normalized USGS hard-reset source locations."""
    return _normalize_replace_relocated_points(corrections)
