"""Explicit USGS station mappings used by source/sink corrections."""

from pathlib import Path

import yaml

USGS_FLOW_STATION_BY_NAME = {
    # Forced/restored artificial-island sources
    "Wando": None,
    "Turkey": "02172035",
    "Buffalo Bluff": "02244040",
    "Dunns Creek": "02244440",

    # Existing sources with explicit observation overrides
    "Delaware": "01463500",
    "Hudson River": "01358000",
}

USGS_TEMPERATURE_STATION_BY_NAME = {
    "Wando": None,
    "Turkey": None,
    "Buffalo Bluff": "02244040",
    "Dunns Creek": "02244440",
    "Delaware": "01463500",

    # Hudson River at Albany, NY
    "Hudson River": "01359139",
}

USGS_FLOW_PARAMETER_ID = "00060"
USGS_TEMPERATURE_PARAMETER_ID = "00010"
CFS_TO_CMS = 0.028316846592

def _load_feature_station_links(section: str) -> dict[int, str]:
    """Read one scope of manual NWM-to-USGS associations."""
    config_file = Path(__file__).with_name("manual_nwm2usgs.yml")
    with config_file.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)

    links = config[section]
    if not isinstance(links, dict):
        raise ValueError(f"{section} must be a mapping in {config_file}")

    result = {}
    for feature_id, station_id in links.items():
        if not isinstance(feature_id, int) or not isinstance(station_id, str):
            raise ValueError(
                f"{section} requires integer FeatureIDs and quoted "
                f"string station IDs in {config_file}"
            )
        result[feature_id] = station_id
    return result


NWM_TO_USGS_FLOW_ADJUSTMENT = _load_feature_station_links("flow_adjustment")
NWM_TO_USGS_TEMPERATURE_SEARCH = _load_feature_station_links(
    "temperature_search_additions"
)
