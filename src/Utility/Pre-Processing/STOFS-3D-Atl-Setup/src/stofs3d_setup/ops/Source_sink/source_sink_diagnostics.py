"""Diagnostic writers for SCHISM source/sink forcing."""

from pathlib import Path

import numpy as np

from pylib_experimental.schism_file import source_sink
from stofs3d_setup.ops.Source_sink.source_sink_components import _data_array
from stofs3d_setup.ops.Source_sink.spatial_selection import (
    _compute_grid_centers,
)


def _write_diagnostics(
    output_dir: Path,
    patched_ss: source_sink,
    hgrid,
) -> None:
    """Write patched source/sink element coordinates and mean flows."""
    output_dir.mkdir(parents=True, exist_ok=True)
    xctr, yctr = _compute_grid_centers(hgrid)

    source_eles = [
        int(ele) for ele in np.asarray(patched_ss.source_eles).reshape(-1)
    ]
    if patched_ss.vsource is not None and source_eles:
        source_mean = np.mean(_data_array(patched_ss.vsource), axis=0)
        source_xyz = np.c_[
            xctr[np.asarray(source_eles) - 1],
            yctr[np.asarray(source_eles) - 1],
            source_mean,
        ]
        np.savetxt(
            output_dir / "patched_vsource.xyz",
            source_xyz,
            fmt="%.10f %.10f %.8f",
            header="lon lat mean_vsource",
            comments="",
        )

    sink_eles_raw = getattr(patched_ss, "sink_eles", None)
    sink_eles = (
        []
        if sink_eles_raw is None
        else [int(ele) for ele in np.asarray(sink_eles_raw).reshape(-1)]
    )
    if patched_ss.vsink is not None and sink_eles:
        sink_mean = np.mean(_data_array(patched_ss.vsink), axis=0)
        sink_xyz = np.c_[
            xctr[np.asarray(sink_eles) - 1],
            yctr[np.asarray(sink_eles) - 1],
            sink_mean,
        ]
        np.savetxt(
            output_dir / "patched_vsink.xyz",
            sink_xyz,
            fmt="%.10f %.10f %.8f",
            header="lon lat mean_vsink",
            comments="",
        )
