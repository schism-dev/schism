"""Lightweight integrity checks for source/sink workflow stages."""

from pathlib import Path

import numpy as np

from pylib_experimental.schism_file import source_sink


def validate_source_sink_stage(ss, stage):
    """Raise with the stage name when forcing is structurally or numerically invalid."""
    if ss.vsource is None and ss.vsink is None:
        raise ValueError(f"{stage}: source and sink are both absent")

    histories = [("vsource", ss.vsource), ("vsink", ss.vsink)]
    if ss.vsource is not None:
        if ss.msource is None or len(ss.msource) != 2:
            raise ValueError(f"{stage}: expected temperature and salinity tracers")
        histories.extend((f"msource[{i}]", tracer) for i, tracer in enumerate(ss.msource))

    for name, history in histories:
        if history is None:
            continue
        if history.df.empty or not history.df.index.is_monotonic_increasing:
            raise ValueError(f"{stage}: {name} has empty or unordered time records")
        if history.df.index.has_duplicates or history.df.columns.has_duplicates:
            raise ValueError(f"{stage}: {name} has duplicate times or elements")
        if not np.isfinite(history.df.to_numpy(copy=False)).all():
            raise ValueError(f"{stage}: {name} has non-finite values")

    if ss.vsink is not None and (ss.vsink.df.to_numpy(copy=False) > 0).any():
        raise ValueError(f"{stage}: vsink contains positive flow")
    if ss.vsource is None:
        return

    source_elements = ss.vsource.df.columns
    for name, tracer in histories[2:]:
        if tracer is None or not tracer.df.columns.equals(source_elements):
            raise ValueError(f"{stage}: {name} elements do not match vsource")

    if (ss.vsource.df.to_numpy(copy=False) < 0).any():
        raise ValueError(f"{stage}: vsource contains negative flow")

    temperature = ss.msource[0].df.to_numpy(copy=False)
    if ((temperature < -2) & (temperature != -9999)).any():
        raise ValueError(f"{stage}: temperature contains values below -2 other than -9999")
    ambient = (temperature == -9999)
    if (ambient.any(axis=0) & ~ambient.all(axis=0)).any():
        raise ValueError(f"{stage}: temperature mixes -9999 with real values")
    if (ss.msource[1].df.to_numpy(copy=False) < 0).any():
        raise ValueError(f"{stage}: salinity contains negative values")


def validate_source_sink_directory(directory, stage):
    """Load and validate an on-disk text forcing group and its source mapping."""
    directory = Path(directory)
    for name in ("source_sink.in", "sources.json"):
        if not (directory / name).is_file():
            raise ValueError(f"{stage}: missing {directory / name}")
    ss = source_sink.from_files(source_dir=str(directory))
    validate_source_sink_stage(ss, stage)
    return ss
