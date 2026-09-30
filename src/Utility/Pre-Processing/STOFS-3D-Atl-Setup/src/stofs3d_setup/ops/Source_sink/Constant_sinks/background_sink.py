"""Post-processing operations for generated background sinks."""

from pathlib import Path

import numpy as np

from pylib_experimental.schism_file import source_sink
from stofs3d_setup.ops.Source_sink.source_sink_components import (
    _build_timehistory,
    _data_array,
    _time_array,
)


def remove_overlapping_background_sinks(
    background_ss: source_sink,
    base_ss: source_sink,
    output_dir: str | Path | None = None,
) -> source_sink:
    """
    Remove background-sink columns whose elements are already used by base_ss.

    This preserves the original set_constant_sink() function and applies the
    overlap removal afterward.

    Both base source elements and base sink elements are excluded from the
    background sink object.
    """
    if background_ss is None:
        raise ValueError("background_ss must not be None")
    if base_ss is None:
        raise ValueError("base_ss must not be None")

    base_source_eles = {
        int(ele) for ele in np.asarray(base_ss.source_eles).reshape(-1)
    }

    base_sink_raw = getattr(base_ss, "sink_eles", None)
    base_sink_eles = (
        set()
        if base_sink_raw is None
        else {int(ele) for ele in np.asarray(base_sink_raw).reshape(-1)}
    )

    excluded = base_source_eles | base_sink_eles

    background_sink_raw = getattr(background_ss, "sink_eles", None)
    background_sink_eles = (
        []
        if background_sink_raw is None
        else [
            int(ele)
            for ele in np.asarray(background_sink_raw).reshape(-1)
        ]
    )

    if background_ss.vsink is None or not background_sink_eles:
        return background_ss

    sink_time = _time_array(background_ss.vsink)
    sink_data = _data_array(background_ss.vsink)

    keep = np.array(
        [ele not in excluded for ele in background_sink_eles],
        dtype=bool,
    )

    kept_eles = [
        ele
        for ele, keep_it in zip(background_sink_eles, keep)
        if keep_it
    ]
    kept_data = sink_data[:, keep]

    removed = [
        ele
        for ele, keep_it in zip(background_sink_eles, keep)
        if not keep_it
    ]

    print(
        "[BACKGROUND SINK PATCH] removed "
        f"{len(removed)} overlapping background sink element(s)."
    )

    patched_background_ss = source_sink(
        vsource=None,
        vsink=_build_timehistory(sink_time, kept_data, kept_eles),
        msource=None,
    )

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        patched_background_ss.writer(str(output_dir))

    return patched_background_ss
