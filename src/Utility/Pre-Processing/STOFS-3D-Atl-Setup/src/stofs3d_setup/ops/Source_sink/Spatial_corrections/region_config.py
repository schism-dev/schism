"""Normalize configured SCHISM source-zeroing regions."""

from pathlib import Path

from stofs3d_setup.ops.Source_sink.correction_config import _as_dict


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


def zero_source_regions(
    corrections: dict,
    config_dir: str | Path,
) -> list[dict]:
    """Return normalized regions used by the source-zeroing stage."""
    return _normalize_zero_source_regions(
        corrections,
        yaml_dir=Path(config_dir),
    )
