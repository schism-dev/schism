"""Physics core for precipitation partitioning (pure pre‑processor).

Implements a simplified Noah‑MP–style bucket per SCHISM element to split
liquid precipitation into (1) surface runoff and (2) percolation to the
unconfined bucket. Infiltration is limited by both storage and rate caps,
scaled by impervious fraction and a freeze gate.

This module is **framework‑agnostic**: all functions operate on numpy arrays
or xarray DataArrays with a shared broadcasting shape, typically (time, elem).

Key equations (per the project summary)
--------------------------------------
P_liq          = QRAIN + d(ACSNOM)/dt
S              = sum_k (smcmax_k - SOIL_M_k) * DZS_k
I_storage_cap  = S / Δt
I_rate_cap     = dksat_top
q_perc,UG      = d(UGDRNOFF)/dt
I_cap          = (1 - imperv) * f_freeze * min(I_rate_cap, I_storage_cap + q_perc,UG)
I              = min(P_liq, I_cap)
q_surf         = P_liq - I
q_perc         = q_perc,UG

All fluxes in m s^-1; depths in meters; θ are volumetric fractions [0–1].

Notes
-----
* `q_perc` is taken directly from the NWM LDASOUT tendency (UGDRNOFF),
  ensuring mass conservation with the land model. If UGDRNOFF is provided
  as an accumulated depth [m], use `finite_diff_accum()` first.
* Freeze gating can be binary (`freeze_gate_binary`) or smooth
  (`freeze_gate_smooth`). The smooth gate is recommended (α≈3).
* Storage uses the *control depth* defined by DZS per soil layer; typically
  the top ~0.4 m (e.g., [0.1, 0.3]).

"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple
import numpy as np

try:  # optional dependency, but plays well with xarray
    import xarray as xr  # type: ignore
except Exception:  # pragma: no cover
    xr = None  # type: ignore


# -----------------------------------------------------------------------------
# Data structures
# -----------------------------------------------------------------------------

@dataclass
class SoilParams:
    """Per‑element soil/static parameters for the control depth.

    Attributes
    ----------
    smcmax : array_like
        Porosity (θ_sat) for each layer in control depth; shape (layer, elem) or (elem,) if single layer.
    DZS : array_like
        Layer thicknesses [m] participating in storage; shape (layer,) or (layer, elem) for spatially varying thickness.
    dksat_top : array_like
        Saturated hydraulic conductivity [m s^-1] for the *top* layer; rate cap.
    imperv : array_like
        Impervious fraction [0–1] per element.
    """

    smcmax: np.ndarray
    DZS: np.ndarray
    dksat_top: np.ndarray
    imperv: np.ndarray


# -----------------------------------------------------------------------------
# Helper utilities
# -----------------------------------------------------------------------------

def as_array(a):
    """Convert input to numpy array without copying when possible."""
    if xr is not None and isinstance(a, xr.DataArray):
        return a.data
    return np.asarray(a)


def finite_diff_accum(accum: np.ndarray, dt: np.ndarray | float) -> np.ndarray:
    """Convert an accumulated depth time series [m] to a rate [m s^-1].

    Handles hourly/irregular cadences (dt can be scalar or time‑varying).
    Resets (e.g., NWM daily zeroing) are supported by non‑negative differencing.

    Parameters
    ----------
    accum : (time, ...) array
        Accumulated depth increasing with time, occasionally resetting to ~0.
    dt : float or (time,) array
        Time step(s) in seconds between samples.

    Returns
    -------
    rate : (time, ...) array
        Instantaneous rate aligned to the input timestamps; the first step uses forward diff.
    """
    A = as_array(accum)
    # ensure time is axis 0
    dA = np.diff(A, axis=0, prepend=A[[0], ...])
    # guard against resets / negative diffs
    dA = np.maximum(dA, 0.0)

    if np.isscalar(dt):
        return dA / float(dt)
    dt_arr = as_array(dt).reshape((-1,) + (1,) * (A.ndim - 1))
    return dA / dt_arr


# -----------------------------------------------------------------------------
# Freeze gating
# -----------------------------------------------------------------------------

def freeze_gate_binary(SNEQV, SOIL_M, SOIL_W, gate_when_ice: float = 0.1) -> np.ndarray:
    """Binary freeze gate per project summary.

    gate = 0.1 if (SNEQV > 0 or SOIL_M - SOIL_W > 0.02) else 1.0
    """
    sne = as_array(SNEQV)
    sm = as_array(SOIL_M)
    sw = as_array(SOIL_W)
    iced = (sne > 0) | ((sm - sw) > 0.02)
    return np.where(iced, gate_when_ice, 1.0)


def freeze_gate_smooth(SOIL_M, SOIL_W, alpha: float = 3.0, eps: float = 1e-6) -> np.ndarray:
    """Smooth freeze gate based on ice fraction proxy.

    f_freeze = max( (SOIL_M - SOIL_W) / (SOIL_M + eps), 0 )
    gate = (1 - f_freeze) ** alpha
    Bounded in [0, 1].
    """
    sm = as_array(SOIL_M)
    sw = as_array(SOIL_W)
    frac = np.maximum((sm - sw) / (np.abs(sm) + eps), 0.0)
    gate = (1.0 - frac) ** alpha
    return np.clip(gate, 0.0, 1.0)


# -----------------------------------------------------------------------------
# Caps and partitioning
# -----------------------------------------------------------------------------

def storage_capacity(SOIL_M_layers: np.ndarray, soil: SoilParams) -> np.ndarray:
    """Compute storage depth [m] available in control depth.

    Parameters
    ----------
    SOIL_M_layers : array, shape (layer, time, elem) or (layer, elem)
        Volumetric water content per layer within the control depth.
    soil : SoilParams
        Soil/static parameters including smcmax and DZS.

    Returns
    -------
    S : array broadcastable to (time, elem)
        Storage depth available to fill before hitting porosity.
    """
    smcmax = as_array(soil.smcmax)
    DZS = as_array(soil.DZS)
    theta = as_array(SOIL_M_layers)

    # Align shapes: we want (layer, time, elem)
    if theta.ndim == 2:  # (layer, elem) static in time
        theta = theta[:, None, :]
    if smcmax.ndim == 1:
        smcmax = smcmax[:, None]
    if DZS.ndim == 1:
        DZS = DZS[:, None]

    # (layer, time, elem)
    # expand smcmax to match theta shape
    smcmax = smcmax.reshape(len(DZS[0]), 1, -1)
    smcmax = np.repeat(smcmax, theta.shape[1], axis=1)
    deficit = np.maximum(smcmax - theta, 0.0)
    S = np.einsum('lte,l->te', deficit, DZS[0])  # assuming Dz is homogeneous in elem and time
    return S


def infiltration_caps(
    dt: np.ndarray | float,
    SOIL_M_layers: np.ndarray,
    soil: SoilParams,
    qperc_ug: np.ndarray,
    freeze_gate: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute (I_rate_cap, I_storage_cap_effective) as arrays [m s^-1].

    I_rate_cap = dksat_top
    I_storage_cap_effective = (S/Δt + qperc_ug) * (1 - imperv) * freeze_gate
    The final I_cap = min(I_rate_cap, I_storage_cap_effective).
    """
    dt_sec = float(dt) if np.isscalar(dt) else as_array(dt).reshape((-1, 1))

    S = storage_capacity(SOIL_M_layers, soil)  # (time, elem)
    I_storage_cap = S / dt_sec + as_array(qperc_ug)

    scale = (1.0 - as_array(soil.imperv)) * as_array(freeze_gate)
    # Broadcast to (time, elem)
    while scale.ndim < I_storage_cap.ndim:
        scale = np.expand_dims(scale, 0)
    I_storage_cap_eff = I_storage_cap * scale

    I_rate_cap = as_array(soil.dksat_top)
    # Broadcast to (time, elem)
    while I_rate_cap.ndim < I_storage_cap_eff.ndim:
        I_rate_cap = np.expand_dims(I_rate_cap, 0)

    return I_rate_cap, I_storage_cap_eff


def partition(
    P_liq: np.ndarray,
    I_rate_cap: np.ndarray,
    I_storage_cap_eff: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Partition liquid input into infiltration, surface runoff, and percolation.

    Parameters
    ----------
    P_liq : (time, elem)
        Liquid water input rate [m s^-1] = QRAIN + d(ACSNOM)/dt
    I_rate_cap : (time, elem)
        Rate limitation [m s^-1]
    I_storage_cap_eff : (time, elem)
        Storage+percolation limited cap after scaling [m s^-1]

    Returns
    -------
    I_rate : (time, elem)
        Infiltration rate [m s^-1]
    q_surf : (time, elem)
        Surface runoff [m s^-1]
    I_cap : (time, elem)
        The applied infiltration capacity (min of caps) for diagnostics.
    """
    I_cap = np.minimum(I_rate_cap, I_storage_cap_eff)
    I_rate = np.minimum(as_array(P_liq), I_cap)
    q_surf = as_array(P_liq) - I_rate
    return I_rate, q_surf, I_cap


# -----------------------------------------------------------------------------
# High‑level driver
# -----------------------------------------------------------------------------

def compute_fluxes(
    dt: np.ndarray | float,
    QRAIN: np.ndarray,
    ACSNOM_accum: np.ndarray | None,
    UGDRNOFF_accum: np.ndarray | None,
    SOIL_M_layers: np.ndarray,
    SOIL_W_layers: np.ndarray | None,
    SNEQV: np.ndarray | None,
    soil: SoilParams,
    use_smooth_freeze: bool = True,
    alpha: float = 3.0,
    eps: float = 1e-6,
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """End‑to‑end flux computation for one control depth.

    Parameters
    ----------
    dt : float or (time,) array
        Time step(s) in seconds.
    QRAIN : (time, elem)
        Liquid precipitation rate [m s^-1].
    ACSNOM_accum : (time, elem) or None
        Accumulated snowmelt depth [m]. If None, treated as zeros.
    UGDRNOFF_accum : (time, elem) or None
        Accumulated unconfined groundwater runoff depth [m]. If None, q_perc=0.
    SOIL_M_layers : (layer, time, elem)
        Volumetric soil moisture per layer [m^3 m^-3].
    SOIL_W_layers : (layer, time, elem) or None
        Liquid water content per layer (excludes ice). If None, assumes no ice.
    SNEQV : (time, elem) or None
        Snow water equivalent [m]. If None, assumes 0.
    soil : SoilParams
        Static soil parameters (porosity, thicknesses, ks, imperv).

    Returns
    -------
    q_surf : (time, elem)
        Surface runoff [m s^-1]
    q_perc : (time, elem)
        Percolation (from UGDRNOFF) [m s^-1]
    diag : dict
        Diagnostics including P_liq, I, I_cap, freeze_gate, I_rate_cap, I_storage_cap_eff.
    """
    dt_arr = dt

    # P_liq = QRAIN + d(ACSNOM)/dt
    if ACSNOM_accum is None:
        dACSNOM_dt = 0.0
    else:
        dACSNOM_dt = finite_diff_accum(ACSNOM_accum, dt_arr)
    P_liq = as_array(QRAIN) + as_array(dACSNOM_dt)

    # q_perc from UGDRNOFF
    if UGDRNOFF_accum is None:
        q_perc = np.zeros_like(P_liq)
    else:
        q_perc = finite_diff_accum(UGDRNOFF_accum, dt_arr)

    # Freeze gate
    if use_smooth_freeze:
        if SOIL_W_layers is None:
            # No ice info → no gating
            gate = np.ones_like(P_liq)
        else:
            gate = freeze_gate_smooth(SOIL_M_layers[0], SOIL_W_layers[0], alpha=alpha, eps=eps)
    else:
        sne = np.zeros_like(P_liq) if SNEQV is None else SNEQV
        if SOIL_W_layers is None:
            # fallback: treat SOIL_W≈SOIL_M
            gate = freeze_gate_binary(sne, SOIL_M_layers[0], SOIL_M_layers[0])
        else:
            gate = freeze_gate_binary(sne, SOIL_M_layers[0], SOIL_W_layers[0])

    I_rate_cap, I_storage_cap_eff = infiltration_caps(
        dt=dt_arr,
        SOIL_M_layers=SOIL_M_layers,
        soil=soil,
        qperc_ug=q_perc,
        freeze_gate=gate,
    )

    I, q_surf, I_cap = partition(P_liq, I_rate_cap, I_storage_cap_eff)

    diag = {
        "P_liq": P_liq,
        "I": I,
        "I_cap": I_cap,
        "freeze_gate": gate,
        "I_rate_cap": I_rate_cap,
        "I_storage_cap_eff": I_storage_cap_eff,
    }

    return q_surf, q_perc, diag


# -----------------------------------------------------------------------------
# Minimal self‑test (dev aid)
# -----------------------------------------------------------------------------
if __name__ == "__main__":  # simple sanity check
    nt, ne, nl = 5, 3, 2
    dt = 3600.0
    QRAIN = np.full((nt, ne), 2e-7)  # ~0.72 mm/hr
    ACS = np.zeros((nt, ne))
    UG = np.linspace(0, 5e-3, nt).reshape(nt, 1) * np.ones((1, ne))  # accum 0..5 mm
    SOIL_M = np.stack([
        np.full((nt, ne), 0.20),
        np.full((nt, ne), 0.25),
    ], axis=0)
    SOIL_W = SOIL_M.copy()

    soil = SoilParams(
        smcmax=np.array([0.45, 0.45]),
        DZS=np.array([0.1, 0.3]),
        dksat_top=np.full(ne, 1e-6),
        imperv=np.zeros(ne),
    )

    qs, qp, dg = compute_fluxes(dt, QRAIN, ACS, UG, SOIL_M, SOIL_W, None, soil)
    print("q_surf mean (mm/hr):", np.mean(qs) * 3600 * 1000)
    print("q_perc mean (mm/hr):", np.mean(qp) * 3600 * 1000)
    print("I mean (mm/hr):", np.mean(dg["I"]) * 3600 * 1000)
