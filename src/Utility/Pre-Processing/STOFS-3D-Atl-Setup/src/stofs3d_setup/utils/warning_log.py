"""Write user-actionable warnings to the setup run's master warning file."""

from collections.abc import Iterable
from pathlib import Path


def initialize_warning_log(path: str | Path) -> Path:
    """Create an empty warning log for a new formal input-generation run."""
    log_path = Path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        "Warnings recorded during STOFS-3D-ATL input generation:\n",
        encoding="utf-8",
    )
    return log_path


def append_input_warning(
    path: str | Path | None,
    *,
    message: str,
    input_files: Iterable[str],
) -> None:
    """Append one warning and the generated inputs it may affect."""
    if path is None:
        return

    files = tuple(input_files)
    with Path(path).open("a", encoding="utf-8") as log:
        for input_file in files:
            log.write(
                f"WARNING: {message} [generated input: {input_file}]\n"
            )
