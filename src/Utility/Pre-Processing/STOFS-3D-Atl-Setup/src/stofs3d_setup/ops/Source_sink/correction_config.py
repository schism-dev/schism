"""Loading and normalization for configured source/sink corrections."""

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    USGS_FLOW_STATION_BY_NAME,
    USGS_TEMPERATURE_STATION_BY_NAME,
)


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


def _normalize_zero_source_regions(
    patch_info: dict,
    yaml_dir: Path,
) -> list[dict]:
    """Normalize SCHISM regions whose existing source flows are zeroed."""
    raw_regions = patch_info.get("zero_source_regions") or []
    if isinstance(raw_regions, (str, Path, dict)):
        raw_regions = [raw_regions]

    regions = []
    for raw in raw_regions:
        if isinstance(raw, (str, Path)):
            entry = {
                "name": Path(raw).stem,
                "region_file": raw,
            }
        else:
            entry = _as_dict(raw)

        region_value = entry.get(
            "region_file",
            entry.get("rgn_file", entry.get("file")),
        )
        if region_value is None:
            raise ValueError(
                "Each zero_source_regions entry requires region_file: "
                f"{entry}"
            )

        region_file = Path(region_value).expanduser()
        if not region_file.is_absolute():
            region_file = yaml_dir / region_file
        region_file = region_file.resolve()

        if not region_file.is_file():
            raise FileNotFoundError(
                "Zero-source SCHISM region does not exist: "
                f"{region_file}"
            )

        regions.append(
            {
                **entry,
                "name": str(entry.get("name", region_file.stem)),
                "region_file": region_file,
            }
        )

    return regions


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
