"""Explicit USGS station mappings used by source/sink corrections."""

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

# Manual NWM feature-to-USGS associations used by station discovery.
MANUAL_NWM_TO_USGS_FLOW = {
    19406836: "07381490",
    15708755: "02489500",
    18928090: "07375175",
    19269176: "07374000",
    16665157: "02244040",
    2590217: "01463500",
    6186156: "01358000",
}
