"""Shared loading helpers for stage-specific source/sink YAML files."""

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


def _as_dict(value: Any) -> dict:
    """Convert a Pydantic model or mapping-like value to a plain dictionary."""
    if hasattr(value, "model_dump"):
        return value.model_dump()
    return dict(value)


def load_source_sink_corrections(config_file: str | Path | dict) -> dict:
    """Read one stage's settings from YAML or copy an already loaded mapping."""
    if isinstance(config_file, dict):
        return deepcopy(config_file)

    config_file = Path(config_file)
    if not config_file.exists():
        raise FileNotFoundError(
            f"Source/sink correction YAML does not exist: {config_file}"
        )

    with config_file.open("r", encoding="utf-8") as stream:
        data = yaml.safe_load(stream) or {}

    if not isinstance(data, dict):
        raise ValueError(
            "Source/sink correction YAML must contain a top-level mapping: "
            f"{config_file}"
        )
    return data
