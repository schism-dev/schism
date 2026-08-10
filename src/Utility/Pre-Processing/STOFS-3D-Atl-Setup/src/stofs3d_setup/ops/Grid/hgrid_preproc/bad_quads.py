"""Bad-quad splitting API.

The implementation remains with the improvement algorithms because it shares
their SCHISM grid readers and byte-compatible output formatting.
"""

from .improve import fix_bad_quads

__all__ = ["fix_bad_quads"]
