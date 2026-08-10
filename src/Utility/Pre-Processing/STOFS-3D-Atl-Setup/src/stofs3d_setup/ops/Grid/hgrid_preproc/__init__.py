"""STOFS-3D horizontal-grid preprocessing workflow."""

from .constants import STOFS_BOUNDARIES
from .partition_check import dummy_vgrid_text, prepare_partition_check

__all__ = [
    "STOFS_BOUNDARIES",
    "dummy_vgrid_text",
    "fix_bad_quads",
    "prepare_partition_check",
]


def __getattr__(name):
    """Load the dependency-heavy mesh implementation only when requested."""
    if name == "fix_bad_quads":
        from .bad_quads import fix_bad_quads

        return fix_bad_quads
    raise AttributeError(name)
