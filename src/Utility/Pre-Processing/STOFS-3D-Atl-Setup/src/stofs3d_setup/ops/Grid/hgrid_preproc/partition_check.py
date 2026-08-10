"""Prepare—but do not execute—the offline SCHISM partition check."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

from .constants import (
    DUMMY_VGRID_H_S,
    DUMMY_VGRID_LEVELS,
    DUMMY_VGRID_MAX_DEPTH,
    DUMMY_VGRID_THETA_B,
    DUMMY_VGRID_THETA_F,
)


def dummy_vgrid_text() -> str:
    """Return the standard 70-level dummy SZ-coordinate vgrid."""
    lines = [
        "2\n",
        f"{DUMMY_VGRID_LEVELS} 1 {DUMMY_VGRID_MAX_DEPTH:.1f}\n",
        "Z levels \n",
        f"1 {-DUMMY_VGRID_MAX_DEPTH:.1f} \n",
        "S levels \n",
        f"{DUMMY_VGRID_H_S:.1f} {DUMMY_VGRID_THETA_B:.2f} "
        f"{DUMMY_VGRID_THETA_F:.1f}\n",
    ]
    denominator = DUMMY_VGRID_LEVELS - 1
    for level in range(1, DUMMY_VGRID_LEVELS + 1):
        sigma = (level - DUMMY_VGRID_LEVELS) / denominator
        lines.append(f"{level} {sigma:.9f}\n")
    return "".join(lines)


def write_dummy_vgrid(path: Path, overwrite: bool = False) -> None:
    """Atomically create the dummy ``vgrid.in`` without copying a template."""
    if (path.exists() or path.is_symlink()) and not overwrite:
        raise FileExistsError(f"refusing to replace existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", newline="") as stream:
            stream.write(dummy_vgrid_text())
        temporary_path.chmod(0o644)
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)


def _replace_symlink(link: Path, target: Path, overwrite: bool = False) -> None:
    if link.exists() or link.is_symlink():
        if not overwrite:
            raise FileExistsError(f"refusing to replace existing path: {link}")
        link.unlink()
    link.parent.mkdir(parents=True, exist_ok=True)
    relative_target = os.path.relpath(target, start=link.parent)
    link.symlink_to(relative_target)


def prepare_partition_check(
    input_grid: Path,
    output_dir: Path,
    partition_script: Path = Path("~/bin/partition_offline.pl"),
    partition_count: int = 3994,
    overwrite: bool = False,
) -> None:
    """Create fixed-name partition inputs and instructions without running it."""
    hgrid_link = output_dir / "hgrid.gr3"
    vgrid_path = output_dir / "vgrid.in"
    instructions = output_dir / "RUN_PARTITION.md"
    existing = [
        path for path in (hgrid_link, vgrid_path, instructions)
        if path.exists() or path.is_symlink()
    ]
    if existing and not overwrite:
        joined = "\n  ".join(str(path) for path in existing)
        raise FileExistsError(
            "partition-check output already exists; enable overwrite to replace it:\n  "
            + joined
        )
    if partition_count < 1:
        raise ValueError("partition_count must be positive")

    output_dir.mkdir(parents=True, exist_ok=True)
    _replace_symlink(hgrid_link, input_grid, overwrite=overwrite)
    write_dummy_vgrid(vgrid_path, overwrite=overwrite)

    command = f"{partition_script.expanduser().absolute()} {partition_count}"
    instructions.write_text(
        "# Standalone offline partition check\n\n"
        "The preprocessing orchestrator deliberately does not run this "
        "memory-intensive step.\n"
        "From this directory, run:\n\n"
        "```bash\n"
        f"cd {output_dir}\n"
        f"{command}\n"
        "```\n\n"
        "Inputs:\n\n"
        f"- `hgrid.gr3` -> `{os.readlink(hgrid_link)}`\n"
        "- `vgrid.in`: generated 70-level dummy SZ grid\n",
        encoding="utf-8",
    )
    print(f"Generated dummy vgrid: {vgrid_path}")
    print(f"Partition input link: {hgrid_link} -> {os.readlink(hgrid_link)}")
    print("Partitioning was not run. Execute separately:")
    print(f"  cd {output_dir}")
    print(f"  {command}")
