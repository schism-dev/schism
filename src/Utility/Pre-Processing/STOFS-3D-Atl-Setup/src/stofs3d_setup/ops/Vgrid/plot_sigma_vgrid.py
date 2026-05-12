
from __future__ import annotations
from pathlib import Path
from typing import Union, Iterable, Optional


import numpy as np
import matplotlib.pyplot as plt

def read_vgrid(vgrid_file):
    with open(vgrid_file, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    idx = 0
    ivcor = int(lines[idx]); idx += 1
    if ivcor != 2:
        raise ValueError("Only ivcor=2 supported")

    nvrt, kz, h_s = map(float, lines[idx].split())
    nvrt = int(nvrt)
    kz = int(kz)
    idx += 1

    assert lines[idx].startswith("Z")
    idx += 1
    z_levels = []
    for _ in range(kz):
        _, z = lines[idx].split()
        z_levels.append(float(z))
        idx += 1

    assert lines[idx].startswith("S")
    idx += 1
    theta_f, theta_b, hc = map(float, lines[idx].split())
    idx += 1

    sigma = []
    for _ in range(nvrt - kz + 1):
        _, s = lines[idx].split()
        sigma.append(float(s))
        idx += 1

    return {
        "nvrt": nvrt,
        "kz": kz,
        "z_levels": np.array(z_levels),
        "sigma": np.array(sigma),
        "theta_f": theta_f,
        "theta_b": theta_b,
        "hc": hc,
    }


def compute_vertical_coords(vgrid, depth):
    h = depth
    etal = 0.0

    sigma = vgrid["sigma"]
    theta_f = vgrid["theta_f"]
    theta_b = vgrid["theta_b"]
    hc = vgrid["hc"]

    sinh_tf = np.sinh(theta_f) if theta_f != 0 else 1.0
    tanh_half = np.tanh(theta_f * 0.5)

    cs = (1 - theta_b) * np.sinh(theta_f * sigma) / sinh_tf \
       + theta_b * (np.tanh(theta_f * (sigma + 0.5)) - tanh_half) / (2 * tanh_half)

    z_sigma = etal * (1 + sigma) + hc * sigma + (h-hc) * cs
    return z_sigma


def plot_vgrid(vgrid_file, depths=(10, 50, 200), label_layers=False):
    vgrid = read_vgrid(vgrid_file)

    plt.figure(figsize=(6, 6))

    for d in depths:
        z = compute_vertical_coords(vgrid, d)
        x = np.full_like(z, d, dtype=float)

        # vertical grid line
        plt.plot(x, z, "k-", lw=1)

        # layer positions
        plt.scatter(x, z, s=20, c="r", zorder=3)

        if label_layers:
            for k, zk in enumerate(z, start=1):
                plt.text(d + 0.5, zk, f"{k}", fontsize=8, va="center")

    plt.gca().invert_yaxis()
    plt.xlabel("Water depth (m)")
    plt.ylabel("z (m)")
    plt.title("SCHISM vertical grid (layer positions marked)")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()


def write_schism_vgrid_in(
    nvrt: int,
    outpath: Union[str, Path] = "vgrid.in",
    *,
    # keep template constants fixed per your request
    kz: int = 1,
    hs: float = 5000.0,
    z_level: float = -5000.0,
    theta_f: float = 10.0,
    theta_b: float = 0.0,
    hc: float = 5.0,
    float_fmt: str = "{:.9f}",
    bottom_cluster_p: Optional[float] = None,
) -> Path:
    """
    Generate a SCHISM vgrid.in (ivcor=2) using your template:
      - always ivcor=2
      - always kz=1 and one Z level at -5000 (configurable via args but default fixed)
      - only sigma coordinates vary with nvrt (or user-provided sigma)

    Parameters
    ----------
    nvrt : int
        Total number oGf vertical levels (must be >= 2 since kz=1 and need at least 1 sigma level).
    outpath : str|Path
        Output filename.
    theta_f (float):
        Vertical stretching strength parameter for sigma (S) coordinates.
        Controls the overall nonlinearity of the vertical coordinate
        transformation. Larger values increase resolution near both the
        surface and the bottom, while smaller values approach uniform
        sigma spacing. Typical values range from ~5 (mild stretching)
        to ~12 (strong stretching). Setting theta_f = 0 disables stretching
        and yields a linear sigma grid.

    theta_b (float):
        Bottom intensification parameter for sigma (S) coordinates.
        Redistributes the vertical stretching introduced by theta_f toward
        the bottom boundary. Values range from 0 to 1:
        - theta_b = 0   → symmetric stretching (no bottom preference)
        - theta_b ≈ 0.3–0.6 → moderate bottom refinement
        - theta_b ≈ 0.7–1.0 → strong bottom refinement
        Increasing theta_b increases vertical resolution near the seabed
        at the expense of mid-water resolution.
    
    bottom_cluster_p : float, optional
        If provided and sigma is None, generate bottom-clustered sigma using:
          sigma(s) = -(1 - (1 - s)^p),  s in [0,1]
            p < 1: more clustering near bottom (sigma=-1)
            p = 1: uniform sigma
            p > 1: more clustering near surface (sigma=0)
        If None, default uniform sigma.

    Returns
    -------
    Path to written file.
    """
    if nvrt < 2:
        raise ValueError("nvrt must be >= 2 (kz=1 leaves at least 1 sigma level).")
    if kz != 1:
        raise ValueError("This helper is designed for kz=1 as requested.")
    nsig = nvrt - kz + 1
    if nsig < 1:
        raise ValueError("nvrt too small: must satisfy nvrt - kz >= 1.")
    
    # Linear sigma: i=1 -> -1, i=nvrt -> 0
    if bottom_cluster_p is not None:
        s = np.linspace(0, 1, nsig)
        sigma = -1.0 + (1.0 - (1.0 - s) ** bottom_cluster_p)
    else:
        sigma = -1.0 + (np.arange(nvrt) / (nvrt - 1.0))  # length nvrt

    with open(outpath, "w") as f:
        f.write("2\n")
        f.write(f"{nvrt:d} {kz} {hs}\n")  # total levels, kz=1, h_s=5000
        f.write("Z levels \n")
        f.write(f"{kz} {z_level} \n")
        f.write("S levels \n")
        f.write(f"{hc} {theta_b} {theta_f}\n")
        for i, s in enumerate(sigma, start=1):
            f.write(f"{i} {float_fmt.format(s)}\n")

    return outpath


# Example usage:
# write_schism_vgrid_in(25, "vgrid_25.in")                 # uniform sigma like your template
# write_schism_vgrid_in(41, "vgrid_41.in")                 # uniform sigma for nvrt=41
# write_schism_vgrid_in(41, "vgrid_41_bot.in", bottom_cluster_p=3.0)  # bottom-refined sigma



# ===== run =====
vgrid_file = "/sciclone/schism10/feiye/STOFS3D-v8/I35/vgrid.in"
write_schism_vgrid_in(
    nvrt=70,
    outpath=vgrid_file ,
    kz=1,
    hs=5000.0,
    z_level=-5000.0,
    theta_f=4.0,
    theta_b=0.75,
    hc=5.0,
    bottom_cluster_p=None   #0.8; set to None for uniform sigma
)
plot_vgrid(
    vgrid_file,
    depths=[10, 20, 30, 40],
    label_layers=False   # set True if you want layer indices
)
