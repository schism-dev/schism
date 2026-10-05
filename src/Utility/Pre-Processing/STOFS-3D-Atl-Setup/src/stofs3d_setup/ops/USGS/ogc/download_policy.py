"""Load shared policy for USGS observation downloads."""

from pathlib import Path
from typing import Any

import yaml


def _load_usgs_download_policy() -> tuple[dict[str, Any], float]:
    """Read and validate shared download settings from the adjacent YAML."""
    config_file = Path(__file__).with_name("usgs_download_policy.yml")
    with config_file.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)

    if not isinstance(config, dict) or not isinstance(
        config.get("download"), dict
    ):
        raise ValueError(f"download must be a mapping in {config_file}")

    policy = dict(config["download"])
    primary_days = policy.get("primary_chunk_days")
    retry_days = policy.get("retry_chunk_days")
    cache_prefix = policy.get("cache_prefix")
    max_gap_days = policy.get("max_missing_gap_days")

    if type(primary_days) is not int or primary_days <= 0:
        raise ValueError(
            f"download.primary_chunk_days must be a positive integer in "
            f"{config_file}"
        )
    if not isinstance(retry_days, list) or any(
        type(days) is not int or days <= 0 for days in retry_days
    ):
        raise ValueError(
            f"download.retry_chunk_days must be a list of positive integers "
            f"in {config_file}"
        )
    if any(
        later >= earlier
        for earlier, later in zip([primary_days, *retry_days], retry_days)
    ):
        raise ValueError(
            f"download chunk sizes must be strictly decreasing in {config_file}"
        )
    if not isinstance(cache_prefix, str) or not cache_prefix.strip():
        raise ValueError(
            f"download.cache_prefix must be a non-empty string in "
            f"{config_file}"
        )
    if Path(cache_prefix).name != cache_prefix:
        raise ValueError(
            f"download.cache_prefix must not contain a path in {config_file}"
        )
    if (
        not isinstance(max_gap_days, (int, float))
        or isinstance(max_gap_days, bool)
        or max_gap_days <= 0
    ):
        raise ValueError(
            f"download.max_missing_gap_days must be positive in {config_file}"
        )

    policy.pop("max_missing_gap_days")
    return policy, float(max_gap_days)


USGS_DOWNLOAD_POLICY, USGS_MAX_MISSING_GAP_DAYS = (
    _load_usgs_download_policy()
)
