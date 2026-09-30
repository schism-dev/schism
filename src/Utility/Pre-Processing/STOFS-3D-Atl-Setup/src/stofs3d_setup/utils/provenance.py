"""Reusable provenance records for STOFS-3D preprocessing workflows."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from enum import Enum
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import socket
import subprocess
import sys
from types import TracebackType
from typing import Mapping, Sequence


_COPY_IGNORE_PATTERNS = ("__pycache__", "*.pyc", "*.pyo")
_HASH_CHUNK_BYTES = 1024 * 1024


def _json_value(value: object) -> object:
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if isinstance(value, (date, datetime)):
        return value.isoformat()
    if isinstance(value, Enum):
        return _json_value(value.value)
    if hasattr(value, "model_dump"):
        return _json_value(value.model_dump())
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_value(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(_json_value(value), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(_HASH_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def _file_record(path: Path, hash_file: bool) -> dict[str, object]:
    absolute_path = path.expanduser().absolute()
    record: dict[str, object] = {"path": str(absolute_path)}
    if not absolute_path.exists() and not absolute_path.is_symlink():
        record.update({"exists": False, "type": "missing"})
        return record

    stat_result = (
        absolute_path.lstat() if absolute_path.is_symlink() else absolute_path.stat()
    )
    record.update(
        {
            "exists": True,
            "size_bytes": stat_result.st_size,
            "modified_time_utc": datetime.fromtimestamp(
                stat_result.st_mtime, timezone.utc
            ).isoformat(),
        }
    )
    if absolute_path.is_symlink():
        record["type"] = "symlink"
        record["target"] = os.readlink(absolute_path)
    elif absolute_path.is_dir():
        record["type"] = "directory"
    else:
        record["type"] = "file"
    if hash_file and absolute_path.is_file():
        record["sha256"] = _sha256(absolute_path)
    return record


def _path_manifest(
    paths: Sequence[Path],
    *,
    hash_files: bool,
    excluded_roots: Sequence[Path] = (),
) -> list[dict[str, object]]:
    excluded = [path.expanduser().absolute() for path in excluded_roots]
    records: list[dict[str, object]] = []
    seen: set[Path] = set()

    def is_excluded(path: Path) -> bool:
        return any(path == root or root in path.parents for root in excluded)

    for requested_path in paths:
        absolute_path = requested_path.expanduser().absolute()
        candidates = [absolute_path]
        if absolute_path.is_dir():
            candidates.extend(sorted(absolute_path.rglob("*")))
        for candidate in candidates:
            if candidate in seen or is_excluded(candidate):
                continue
            seen.add(candidate)
            records.append(_file_record(candidate, hash_file=hash_files))
    return records


def _copy_source_paths(
    snapshot_dir: Path,
    source_paths: Sequence[Path],
    source_base: Path,
) -> tuple[list[Path], list[dict[str, object]]]:
    source_base = source_base.expanduser().absolute()
    snapshot_source_root = snapshot_dir / "source"
    relevant_paths: list[Path] = []
    records: list[dict[str, object]] = []

    for source_path in source_paths:
        source = source_path.expanduser().absolute()
        if not source.exists() and not source.is_symlink():
            raise FileNotFoundError(f"provenance source does not exist: {source}")
        try:
            relative_path = source.relative_to(source_base)
        except ValueError as exc:
            raise ValueError(
                f"provenance source {source} is outside source base {source_base}"
            ) from exc
        if relative_path == Path("."):
            raise ValueError("source_base itself cannot be a snapshot source")

        destination = snapshot_source_root / relative_path
        if destination.exists() or destination.is_symlink():
            raise ValueError(f"overlapping provenance sources target {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)

        if source.is_dir():
            shutil.copytree(
                source,
                destination,
                symlinks=True,
                ignore=shutil.ignore_patterns(*_COPY_IGNORE_PATTERNS),
            )
            source_files = sorted(
                path
                for path in source.rglob("*")
                if not any(
                    part == "__pycache__" or part.endswith((".pyc", ".pyo"))
                    for part in path.relative_to(source).parts
                )
                and (path.is_symlink() or not path.is_dir())
            )
        else:
            if source.is_symlink():
                destination.symlink_to(os.readlink(source))
            else:
                shutil.copy2(source, destination)
            source_files = [source]

        for source_file in source_files:
            copied_file = (
                destination / source_file.relative_to(source)
                if source.is_dir()
                else destination
            )
            source_record = _file_record(source_file, hash_file=True)
            copied_record = _file_record(copied_file, hash_file=True)
            if source_record.get("sha256") != copied_record.get("sha256"):
                raise OSError(f"source snapshot verification failed: {source_file}")
            if source_record.get("target") != copied_record.get("target"):
                raise OSError(f"source symlink verification failed: {source_file}")
            records.append(
                {
                    "source": str(source_file),
                    "snapshot": str(copied_file.relative_to(snapshot_dir)),
                    "size_bytes": copied_record.get("size_bytes"),
                    "sha256": copied_record.get("sha256"),
                    "type": copied_record.get("type"),
                    **(
                        {"target": copied_record["target"]}
                        if "target" in copied_record
                        else {}
                    ),
                }
            )
        relevant_paths.append(source)
    return relevant_paths, records


def _git_command(repository: Path, arguments: Sequence[str]) -> str:
    result = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        return f"UNAVAILABLE: {detail}\n"
    return result.stdout


def _write_git_records(snapshot_dir: Path, relevant_paths: Sequence[Path]) -> None:
    search_path = relevant_paths[0]
    if search_path.is_file():
        search_path = search_path.parent
    repository_text = _git_command(search_path, ["rev-parse", "--show-toplevel"])
    if repository_text.startswith("UNAVAILABLE:"):
        (snapshot_dir / "git_info.txt").write_text(
            repository_text, encoding="utf-8"
        )
        (snapshot_dir / "git_diff.patch").write_text(
            repository_text, encoding="utf-8"
        )
        return

    repository = Path(repository_text.strip())
    pathspecs = [str(path.relative_to(repository)) for path in relevant_paths]
    commit = _git_command(repository, ["rev-parse", "HEAD"]).strip()
    branch = _git_command(repository, ["branch", "--show-current"]).strip()
    status = _git_command(repository, ["status", "--short", "--", *pathspecs])
    status_text = status if status else "(clean)\n"
    (snapshot_dir / "git_info.txt").write_text(
        f"repository: {repository}\n"
        f"commit: {commit}\n"
        f"branch: {branch or '(detached HEAD)'}\n"
        "relevant_status:\n"
        f"{status_text}",
        encoding="utf-8",
    )
    diff = _git_command(repository, ["diff", "--binary", "HEAD", "--", *pathspecs])
    (snapshot_dir / "git_diff.patch").write_text(diff, encoding="utf-8")


def _write_environment(path: Path) -> None:
    packages = sorted(
        {
            f"{distribution.metadata.get('Name', 'UNKNOWN')}=={distribution.version}"
            for distribution in importlib.metadata.distributions()
        },
        key=str.casefold,
    )
    lines = [
        f"timestamp_utc: {datetime.now(timezone.utc).isoformat()}",
        f"hostname: {socket.gethostname()}",
        f"platform: {platform.platform()}",
        f"python_executable: {sys.executable}",
        f"python_version: {sys.version.replace(os.linesep, ' ')}",
        "",
        "installed_packages:",
        *(f"  {package}" for package in packages),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_readme(snapshot_dir: Path, workflow_name: str) -> None:
    (snapshot_dir / "README.md").write_text(
        f"# {workflow_name} provenance record\n\n"
        "This directory was created before processing began. It contains the "
        "exact selected source files, resolved parameters, input identities, Git "
        "state, launch command, and Python environment.\n\n"
        "The copied package uses a standard `src/` layout. After recreating the "
        "recorded environment, inspect `command.txt` and run copied modules with "
        "`PYTHONPATH=source/src`.\n\n"
        "`status.json` reports whether the workflow completed, failed through a "
        "handled exception, or stopped before finalization. `outputs.json` is "
        "written only after successful completion.\n",
        encoding="utf-8",
    )


def _create_provenance_record(
    output_root: Path,
    *,
    workflow_name: str,
    parameters: Mapping[str, object],
    command: Sequence[str],
    source_paths: Sequence[Path],
    source_base: Path,
    input_paths: Sequence[Path] = (),
) -> Path:
    """Create and verify a timestamped provenance record before processing."""
    if not source_paths:
        raise ValueError("at least one provenance source path is required")

    output_root = output_root.expanduser().absolute()
    provenance_root = output_root / "Reproduce"
    for source_path in source_paths:
        source = source_path.expanduser().absolute()
        if source.is_dir() and source in provenance_root.parents:
            raise ValueError(
                f"provenance output {provenance_root} cannot be inside "
                f"snapshotted source directory {source}"
            )

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    snapshot_dir = provenance_root / f"run_{timestamp}"
    snapshot_dir.mkdir(parents=True, exist_ok=False)

    relevant_paths, source_records = _copy_source_paths(
        snapshot_dir, source_paths, source_base
    )
    _write_json(snapshot_dir / "source_manifest.json", source_records)
    _write_json(
        snapshot_dir / "inputs.json",
        _path_manifest(input_paths, hash_files=True),
    )
    command_text = f"cd {shlex.quote(str(Path.cwd()))}\n{shlex.join(command)}\n"
    (snapshot_dir / "command.txt").write_text(command_text, encoding="utf-8")
    _write_json(snapshot_dir / "parameters.json", parameters)
    _write_git_records(snapshot_dir, relevant_paths)
    _write_environment(snapshot_dir / "environment.txt")
    _write_readme(snapshot_dir, workflow_name)
    _write_json(
        snapshot_dir / "status.json",
        {
            "status": "started",
            "workflow": workflow_name,
            "started_at_utc": datetime.now(timezone.utc),
        },
    )
    return snapshot_dir


def _complete_provenance_record(
    snapshot_dir: Path,
    output_paths: Sequence[Path],
    *,
    hash_outputs: bool = False,
) -> None:
    """Record successful completion and observable output metadata."""
    _write_json(
        snapshot_dir / "outputs.json",
        _path_manifest(
            output_paths,
            hash_files=hash_outputs,
            excluded_roots=(snapshot_dir.parent,),
        ),
    )
    status = json.loads((snapshot_dir / "status.json").read_text(encoding="utf-8"))
    status.update(
        {
            "status": "completed",
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        }
    )
    _write_json(snapshot_dir / "status.json", status)


def _fail_provenance_record(snapshot_dir: Path, error: BaseException) -> None:
    """Record a handled workflow failure without suppressing the exception."""
    status = json.loads((snapshot_dir / "status.json").read_text(encoding="utf-8"))
    status.update(
        {
            "status": "failed",
            "failed_at_utc": datetime.now(timezone.utc).isoformat(),
            "error_type": type(error).__name__,
            "error_message": str(error),
        }
    )
    _write_json(snapshot_dir / "status.json", status)


@dataclass
class ProvenanceRun:
    """Create, complete, or fail one workflow record as a context manager."""

    output_root: Path
    workflow_name: str
    parameters: Mapping[str, object]
    command: Sequence[str]
    source_paths: Sequence[Path]
    source_base: Path
    input_paths: Sequence[Path] = ()
    output_paths: Sequence[Path] = ()
    hash_outputs: bool = False
    _snapshot_dir: Path | None = field(default=None, init=False, repr=False)

    @property
    def path(self) -> Path:
        """Return the record path after the context has been entered."""
        if self._snapshot_dir is None:
            raise RuntimeError("the provenance context has not been entered")
        return self._snapshot_dir

    def __enter__(self) -> "ProvenanceRun":
        if self._snapshot_dir is not None:
            raise RuntimeError("a ProvenanceRun instance cannot be reused")
        self._snapshot_dir = _create_provenance_record(
            self.output_root,
            workflow_name=self.workflow_name,
            parameters=self.parameters,
            command=self.command,
            source_paths=self.source_paths,
            source_base=self.source_base,
            input_paths=self.input_paths,
        )
        print(f"Provenance record: {self._snapshot_dir}")
        return self

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        del exception_type, traceback
        if exception is not None:
            try:
                _fail_provenance_record(self.path, exception)
            except Exception as provenance_error:
                print(
                    f"WARNING: could not record workflow failure: {provenance_error}",
                    file=sys.stderr,
                )
            return False

        existing_outputs = [
            path
            for path in self.output_paths
            if path.exists() or path.is_symlink()
        ]
        try:
            _complete_provenance_record(
                self.path,
                existing_outputs,
                hash_outputs=self.hash_outputs,
            )
        except BaseException as completion_error:
            try:
                _fail_provenance_record(self.path, completion_error)
            except Exception as provenance_error:
                print(
                    "WARNING: could not record provenance-finalization failure: "
                    f"{provenance_error}",
                    file=sys.stderr,
                )
            raise
        return False

    @classmethod
    def for_hgrid(
        cls,
        args: argparse.Namespace,
        argv: Sequence[str] | None,
    ) -> "ProvenanceRun":
        """Build the standard record for the hgrid preprocessing CLI."""
        project_root = Path(__file__).resolve().parents[3]
        package_root = project_root / "src/stofs3d_setup"
        output_root = Path(args.output_root)
        products = {
            "improve": output_root / "Improve/hgrid.ll",
            "split": output_root / "Split_quads/hgrid.gr3.new",
            "boundary": output_root / "Bnd/hgrid_with_bnd.gr3",
        }
        resume_inputs = {
            "split": products["improve"],
            "boundary": products["split"],
            "partition": products["boundary"],
        }
        input_paths = [Path(args.grid)]
        if args.start_at in resume_inputs:
            input_paths.append(resume_inputs[args.start_at])
        command = (
            list(sys.argv)
            if argv is None
            else [
                sys.executable,
                "-m",
                "stofs3d_setup.ops.Grid.hgrid_preproc.cli",
                *argv,
            ]
        )
        partition_dir = output_root / "Partition_check"
        return cls(
            output_root=output_root,
            workflow_name="STOFS-3D hgrid preprocessing",
            parameters=vars(args),
            command=command,
            source_paths=(
                package_root / "ops/Grid/hgrid_preproc",
                package_root / "__init__.py",
                package_root / "ops/__init__.py",
                package_root / "ops/Grid/__init__.py",
                package_root / "utils/__init__.py",
                package_root / "utils/projection.py",
                package_root / "utils/provenance.py",
                project_root / "pyproject.toml",
            ),
            source_base=project_root,
            input_paths=input_paths,
            output_paths=(
                *products.values(),
                partition_dir / "hgrid.gr3",
                partition_dir / "vgrid.in",
                partition_dir / "RUN_PARTITION.md",
            ),
        )

    @classmethod
    def for_stofs3d_driver(
        cls,
        arguments: Mapping[str, object],
    ) -> "ProvenanceRun":
        """Build the standard record for STOFS-3D Atlantic input generation."""
        required_names = {
            "hgrid_path",
            "vgrid_path",
            "config",
            "project_dir",
            "runid",
            "scr_dir",
            "input_files",
        }
        missing_names = sorted(required_names.difference(arguments))
        if missing_names:
            raise ValueError(
                "missing STOFS-3D driver provenance arguments: "
                + ", ".join(missing_names)
            )

        hgrid_path = Path(os.fspath(arguments["hgrid_path"]))
        project_dir = Path(os.fspath(arguments["project_dir"]))
        runid = str(arguments["runid"])
        vgrid_value = arguments["vgrid_path"]
        project_root = Path(__file__).resolve().parents[3]
        model_input_path = project_dir / f"I{runid}"
        original_configs = sorted(model_input_path.glob("*.yml")) + sorted(
            model_input_path.glob("*.yaml")
        )
        inputs = [hgrid_path, *original_configs]
        if vgrid_value is not None:
            inputs.append(Path(os.fspath(vgrid_value)))
        return cls(
            output_root=model_input_path,
            workflow_name="STOFS-3D Atlantic input generation",
            parameters=dict(arguments),
            command=list(sys.argv),
            source_paths=(
                project_root / "src/stofs3d_setup",
                project_root / "pyproject.toml",
            ),
            source_base=project_root,
            input_paths=inputs,
            output_paths=(model_input_path,),
        )
