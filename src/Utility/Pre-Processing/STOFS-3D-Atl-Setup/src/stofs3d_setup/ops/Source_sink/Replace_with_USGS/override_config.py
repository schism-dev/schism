"""Normalize selected USGS overrides for existing relocated sources."""

from stofs3d_setup.ops.Source_sink.correction_config import _as_dict
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    USGS_FLOW_STATION_BY_NAME,
    USGS_TEMPERATURE_STATION_BY_NAME,
)


def _normalize_replace_relocated_points(patch_info: dict) -> list[dict]:
    """Normalize selected existing-source replacement locations."""
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
        if (
            name not in USGS_FLOW_STATION_BY_NAME
            and name not in USGS_TEMPERATURE_STATION_BY_NAME
        ):
            raise ValueError(
                f"{name!r} is listed for relocated-source replacement, "
                "but no USGS flow or temperature station is configured"
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


def source_override_points(corrections: dict) -> list[dict]:
    """Return normalized explicit USGS source overrides."""
    return _normalize_replace_relocated_points(corrections)
