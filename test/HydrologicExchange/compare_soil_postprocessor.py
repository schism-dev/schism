"""Compare the Fortran soil calculation with the unchanged Python reference."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from io import StringIO

import numpy as np
import xarray as xr


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_PATH = Path(__file__).with_name("soil_parity_fixture.dat")
SOIL_REFERENCE_PATH = (
    REPOSITORY_ROOT
    / "src/Utility/Pre-Processing/STOFS-3D-Atl-Setup/src/stofs3d_setup/ops/Source_sink/Soil"
)
COMPARISON_FIELDS = (
    "q_surface",
    "q_perc",
    "P_liq",
    "I",
    "I_cap",
    "freeze_gate",
    "I_rate_cap",
    "I_storage_cap_eff",
    "available_capacity",
)
ELEMENT_X = np.array([0.25, 0.75, 0.50])
ELEMENT_Y = np.array([0.25, 0.25, 0.75])
ELEMENT_AREA = np.array([100.0, 150.0, 200.0])


def compile_fortran_driver(build_path: Path) -> Path:
    compiler = shlex.split(os.environ.get("FC", "gfortran"))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError("Set FC to an available Fortran compiler, such as gfortran.")

    compiler_name = Path(compiler[0]).name.lower()
    if "gfortran" in compiler_name:
        compiler_flags = ["-cpp", "-ffree-line-length-none"]
        module_flags = ["-J", str(build_path)]
    elif "ifx" in compiler_name or "ifort" in compiler_name:
        compiler_flags = ["-fpp"]
        module_flags = ["-module", str(build_path)]
    else:
        compiler_flags = []
        module_flags = []

    core_source = REPOSITORY_ROOT / "src/Core/schism_glbl.F90"
    soil_source = REPOSITORY_ROOT / "src/HydrologicExchange/hydrologic_exchange_soil.F90"
    driver_source = Path(__file__).with_name("soil_fortran_driver.F90")
    core_object = build_path / "schism_glbl.o"
    executable = build_path / "soil_fortran_driver"

    subprocess.run(
        [
            *compiler,
            *compiler_flags,
            *module_flags,
            "-c",
            str(core_source),
            "-o",
            str(core_object),
        ],
        check=True,
        cwd=build_path,
    )
    subprocess.run(
        [
            *compiler,
            *compiler_flags,
            *module_flags,
            "-I",
            str(build_path),
            str(core_object),
            str(soil_source),
            str(driver_source),
            "-o",
            str(executable),
        ],
        check=True,
        cwd=build_path,
    )
    return executable


def map_to_element_centers(values: np.ndarray) -> np.ndarray:
    field = xr.DataArray(
        values.reshape((-1, 2, 2)),
        dims=("time", "y", "x"),
        coords={"time": np.arange(values.shape[0]), "x": [0.0, 1.0], "y": [0.0, 1.0]},
    )
    mapped = field.interp(
        x=xr.DataArray(ELEMENT_X, dims="nface"),
        y=xr.DataArray(ELEMENT_Y, dims="nface"),
        method="linear",
    )
    return mapped.values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", type=Path, default=FIXTURE_PATH)
    parser.add_argument(
        "--grid-only",
        action="store_true",
        help="Compare NWM-grid fields without the synthetic element mapping check.",
    )
    args = parser.parse_args()

    sys.path.insert(0, str(SOIL_REFERENCE_PATH))
    from physics import SoilParams, compute_fluxes

    with args.fixture.open(encoding="utf-8") as fixture:
        ntime, nelem = map(int, fixture.readline().split())
    records = np.loadtxt(args.fixture, skiprows=1).reshape(ntime, nelem, -1)

    dt = records[:, 0, 2]
    qrain = records[:, :, 3]
    acsnom = records[:, :, 4]
    accet = records[:, :, 5]
    ugdrnoff = records[:, :, 6]
    soil_m = records[:, :, 7][None, :, :]
    soil_w = records[:, :, 8][None, :, :]
    smcmax = records[0, :, 9][None, :]
    dzs = np.array([records[0, 0, 10]])
    dksat_top = records[0, :, 11]
    imperv = records[0, :, 12]
    soil_parameters = SoilParams(
        smcmax=smcmax,
        DZS=dzs,
        dksat_top=dksat_top,
        imperv=imperv,
    )

    q_surface, q_perc, diagnostics = compute_fluxes(
        dt,
        qrain,
        acsnom,
        accet,
        ugdrnoff,
        soil_m,
        soil_w,
        None,
        soil_parameters,
        use_smooth_freeze=True,
    )
    available_capacity = np.maximum(smcmax[0][None, :] - soil_m[0], 0.0) * dzs[0]
    rate_cap = np.broadcast_to(diagnostics["I_rate_cap"], q_surface.shape)
    expected = np.stack(
        [
            q_surface,
            q_perc,
            diagnostics["P_liq"],
            diagnostics["I"],
            diagnostics["I_cap"],
            diagnostics["freeze_gate"],
            rate_cap,
            diagnostics["I_storage_cap_eff"],
            available_capacity,
        ],
        axis=-1,
    )

    with tempfile.TemporaryDirectory(prefix="hydrologic_exchange-pr1-") as temporary_directory:
        build_path = Path(temporary_directory)
        executable = compile_fortran_driver(build_path)
        result = subprocess.run(
            [str(executable), str(args.fixture)],
            check=True,
            capture_output=True,
            text=True,
        )

    rows = np.loadtxt(StringIO(result.stdout))
    actual = rows[:, 2:].reshape(ntime, nelem, len(COMPARISON_FIELDS))
    np.testing.assert_allclose(actual, expected, rtol=1.e-12, atol=1.e-14)

    print(f"Grid-cell comparison passed for {ntime} times and {nelem} cells.")
    if args.grid_only:
        return

    python_mapped = np.stack(
        [map_to_element_centers(q_surface), map_to_element_centers(q_perc)], axis=-1
    )
    fortran_mapped = np.stack(
        [map_to_element_centers(actual[:, :, 0]), map_to_element_centers(actual[:, :, 1])],
        axis=-1,
    )
    np.testing.assert_allclose(fortran_mapped, python_mapped, rtol=1.e-12, atol=1.e-14)

    python_volume = np.sum(map_to_element_centers(q_surface) * ELEMENT_AREA * dt[:, None])
    fortran_volume = np.sum(fortran_mapped[:, :, 0] * ELEMENT_AREA * dt[:, None])
    np.testing.assert_allclose(fortran_volume, python_volume, rtol=1.e-12, atol=1.e-14)

    print("Mapped q_surface and q_perc comparison passed at 3 element centers.")
    print(f"Fixture integrated q_surface volume: {python_volume:.12g} m^3")


if __name__ == "__main__":
    main()
