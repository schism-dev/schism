"""Load USGS station mappings for soft-link source corrections."""

from pathlib import Path
from typing import Any

import yaml


USGS_FLOW_PARAMETER_ID = "00060"
USGS_TEMPERATURE_PARAMETER_ID = "00010"
CFS_TO_CMS = 0.028316846592


def _load_manual_usgs_config() -> tuple[dict[str, Any], Path]:
    """Read the YAML configuration for soft-link USGS corrections."""
    config_file = Path(__file__).with_name("manual_nwm2usgs.yml")
    with config_file.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError(
            f"USGS soft-link configuration must be a mapping: {config_file}"
        )
    return config, config_file


def _load_feature_station_links(
    config: dict[str, Any], section: str, config_file: Path
) -> tuple[dict[int, str], set[str], set[str]]:
    """Load one scope of manual NWM-to-USGS feature associations."""
    links = config.get(section)
    if not isinstance(links, dict):
        raise ValueError(f"{section} must be a mapping in {config_file}")

    result: dict[int, str] = {}
    required: set[str] = set()
    optional: set[str] = set()
    for feature_id, station_config in links.items():
        if type(feature_id) is not int or not isinstance(station_config, dict):
            raise ValueError(
                f"{section} requires integer FeatureIDs and station mappings "
                f"in {config_file}"
            )
        station_id = station_config.get("station_id")
        requirement = station_config.get("download_requirement")
        if not isinstance(station_id, str) or not station_id.strip():
            raise ValueError(
                f"{section}[{feature_id}].station_id must be a quoted "
                f"string in {config_file}"
            )
        if requirement not in {"required", "optional"}:
            raise ValueError(
                f"{section}[{feature_id}].download_requirement must be "
                f"'required' or 'optional' in {config_file}"
            )
        result[feature_id] = station_id
        (required if requirement == "required" else optional).add(station_id)
    return result, required, optional


_CONFIG, _CONFIG_FILE = _load_manual_usgs_config()
TEMPERATURE_SEARCH_DIRECTIONS = ("upstream", "downstream")
(
    NWM_TO_USGS_FLOW_ADJUSTMENT,
    REQUIRED_FLOW_STATIONS,
    OPTIONAL_FLOW_STATIONS,
) = _load_feature_station_links(
    _CONFIG, "flow_adjustment", _CONFIG_FILE
)
(
    NWM_TO_USGS_TEMPERATURE_SEARCH,
    REQUIRED_TEMPERATURE_STATIONS,
    OPTIONAL_TEMPERATURE_STATIONS,
) = _load_feature_station_links(
    _CONFIG, "temperature_search_additions", _CONFIG_FILE
)
