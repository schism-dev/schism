#!/usr/bin/env python3
"""
plot_vqs.py

Diagnostics for VQS vertical grids + transect extraction

What it does:
  1) Read vgrid_master.out and plot the master VQS grids (MATLAB plot_VQS.m equivalent).
  2) If transect.bp exists, generate transect1.out using:
        - hgrid.gr3 (node x,y,dp)
        - vgrid.in  (per-node sigma, ivcor=1)
     by mapping each bp point to the nearest hgrid node (KDTree if available; fallback chunked brute force).
  3) Plot the generated transect1.out.

Files used:
  - vgrid_master.out  (required for master plot)
  - hgrid.gr3         (required if generating transect1.out)
  - vgrid.in          (required if generating transect1.out)
  - transect.bp       (optional; if present, generate transect1.out)

Outputs:
  - transect1.out (if transect.bp exists)
  - vqs_diagnostics.png
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ----------------------------
# hgrid.gr3 reader
# ----------------------------
def read_hgrid_gr3(fname: str | Path):
    """
    Minimal SCHISM hgrid.gr3 reader (mixed tri/quad ok).
    Returns:
      ne, np, x, y, dp, i34, elnode
    """
    fname = Path(fname)
    with fname.open("r") as f:
        _title = f.readline()
        ne, np_ = map(int, f.readline().split()[:2])

        x = np.empty(np_, dtype=float)
        y = np.empty(np_, dtype=float)
        dp = np.empty(np_, dtype=float)

        for i in range(np_):
            parts = f.readline().split()
            x[i] = float(parts[1])
            y[i] = float(parts[2])
            dp[i] = float(parts[3])

        # elements are not needed for transect mapping, but read to be robust
        i34 = np.empty(ne, dtype=int)
        elnode = np.zeros((ne, 4), dtype=int)
        for e in range(ne):
            parts = f.readline().split()
            i34[e] = int(parts[1])
            nodes = list(map(int, parts[2:2 + i34[e]]))
            elnode[e, :i34[e]] = nodes

    return ne, np_, x, y, dp, i34, elnode


# ----------------------------
# vgrid.in reader (ivcor=1)
# ----------------------------
def read_vgrid_in_ivcor1(fname: str | Path):
    """
    Reads SCHISM vgrid.in for ivcor=1 (VQS / sigma-type).
    vgrid.in format (for ivcor=1):
      line1: ivcor
      line2: nvrt
      then np lines:
        i  n_dry  sigma(bottom->surface)  (length = kbp)
      where:
        kbp = nvrt + 1 - n_dry
        sigma list is written as sigma_vqs(kbp:1:-1) i.e. bottom->surface
    Returns:
      nvrt (int)
      sigma_node: list of np arrays (length kbp for node i) in surface->bottom order
      kbp: (np,) int
    """
    fname = Path(fname)
    with fname.open("r") as f:
        ivcor = int(f.readline().split()[0])
        if ivcor != 1:
            raise ValueError(f"Expected ivcor=1 but got {ivcor}")
        nvrt = int(f.readline().split()[0])

        sigma_node = []
        kbp_list = []

        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            # i, n_dry, then sigmas (bottom->surface)
            n_dry = int(parts[1])
            kb = nvrt + 1 - n_dry
            kbp_list.append(kb)

            sig_bot_to_surf = np.array(list(map(float, parts[2:])), dtype=float)
            if sig_bot_to_surf.size != kb:
                # Some writers may wrap lines; in that case you'd need a token-based reader.
                # For typical SCHISM vgrid.in it is one line per node.
                raise ValueError(
                    f"Node line has {sig_bot_to_surf.size} sigmas but kbp={kb}. "
                    f"If your vgrid.in wraps lines, tell me and I'll switch to a robust token reader."
                )

            # convert to surface->bottom order for easier z computation
            sig_surf_to_bot = sig_bot_to_surf[::-1]
            sigma_node.append(sig_surf_to_bot)

    kbp = np.array(kbp_list, dtype=int)
    return nvrt, sigma_node, kbp


def sigma_to_znd(nvrt: int, sigma_node: list[np.ndarray], dp: np.ndarray, eta: float = 0.0):
    """
    Convert sigma (surface->bottom, length kbp) to znd(:,node) with extension beyond bottom.
    z = sigma*(eta+dp) + eta, with eta default 0.
    Returns:
      znd: (nvrt, np) float
    """
    np_ = dp.size
    znd = np.empty((nvrt, np_), dtype=float)

    for i in range(np_):
        sig = sigma_node[i]  # length kbp
        kb = sig.size
        z = sig * (eta + dp[i]) + eta
        znd[:kb, i] = z
        # extend beyond bottom for plotting/output consistency
        znd[kb:, i] = -dp[i]
        # enforce exact bottom if you want
        znd[kb - 1, i] = -dp[i]

    return znd


# ----------------------------
# vgrid_master.out reader
# ----------------------------
def load_vgrid_master(fname: str | Path = "vgrid_master.out"):
    """
    vgrid_master.out row format:
      m  nv(m)  hsm(m)  z_mas(1..nvrt_m)
    Returns:
      m: (m_vqs,)
      nv: (m_vqs,)
      hsm: (m_vqs,)
      z_mas: (nvrt_m, m_vqs)
      zcor_m_masked: (m_vqs, nvrt_m) with NaN beyond nv(m)
    """
    arr = np.loadtxt(fname)
    m = arr[:, 0].astype(int)
    nv = arr[:, 1].astype(int)
    hsm = arr[:, 2].astype(float)
    z_rows = arr[:, 3:]          # (m_vqs, nvrt_m)
    z_mas = z_rows.T             # (nvrt_m, m_vqs)

    zcor_m = z_rows.copy()
    for i in range(zcor_m.shape[0]):
        nm = nv[i]
        if nm < zcor_m.shape[1]:
            zcor_m[i, nm:] = np.nan
    return m, nv, hsm, z_mas, zcor_m


# ----------------------------
# transect.bp handling + nearest node mapping
# ----------------------------
def read_transect_bp(fname: str | Path = "transect.bp"):
    """
    Reads SCHISM bp file:
      line1: comment
      line2: npbp
      then npbp lines: id x y [optional ...]
    Returns:
      xybp: (npbp, 2)
    """
    fname = Path(fname)
    with fname.open("r") as f:
        _ = f.readline()
        npbp = int(f.readline().split()[0])
        xybp = np.empty((npbp, 2), dtype=float)
        for i in range(npbp):
            parts = f.readline().split()
            xybp[i, 0] = float(parts[1])
            xybp[i, 1] = float(parts[2])
    return xybp


def cumulative_distance(xy: np.ndarray):
    """
    xy: (n,2)
    returns s: (n,) cumulative distance with s[0]=0
    """
    dxy = xy[1:] - xy[:-1]
    ds = np.sqrt((dxy ** 2).sum(axis=1))
    s = np.empty(xy.shape[0], dtype=float)
    s[0] = 0.0
    s[1:] = np.cumsum(ds)
    return s


def nearest_nodes(x: np.ndarray, y: np.ndarray, xybp: np.ndarray):
    """
    Finds nearest hgrid node for each bp point.
    Uses scipy.spatial.cKDTree if available; otherwise chunked brute-force.

    Returns:
      imap: (npbp,) int, 0-based node indices
    """
    pts = np.column_stack([x, y])
    q = xybp

    # Try KDTree
    try:
        from scipy.spatial import cKDTree  # type: ignore
        tree = cKDTree(pts)
        dist, idx = tree.query(q, k=1)
        return idx.astype(int)
    except Exception:
        pass

    # Fallback: chunked brute force (memory-safe)
    npbp = q.shape[0]
    imap = np.empty(npbp, dtype=int)
    for i in range(npbp):
        dx = x - q[i, 0]
        dy = y - q[i, 1]
        imap[i] = int(np.argmin(dx * dx + dy * dy))
    return imap


def write_transect1_out(
    fname: str | Path,
    xybp: np.ndarray,
    s: np.ndarray,
    imap: np.ndarray,
    dp: np.ndarray,
    kbp: np.ndarray,
    znd: np.ndarray,
):
    """
    Write transect1.out in the same column order as your Fortran:
      i, kbp(nd), x, y, s(i), dp(nd), znd(:,nd)

    Note: Fortran writes znd(:,nd) length = nvrt_m at that time; here we write nvrt from vgrid.in.
    """
    fname = Path(fname)
    nvrt = znd.shape[0]

    with fname.open("w") as f:
        for i in range(xybp.shape[0]):
            nd = int(imap[i])
            row = [
                f"{i+1:6d}",
                f"{int(kbp[nd]):4d}",
                f"{xybp[i,0]:16.6e}",
                f"{xybp[i,1]:16.6e}",
                f"{s[i]:12.3f}",
                f"{dp[nd]:12.3f}",
            ]
            zvals = " ".join(f"{znd[k, nd]:12.3f}" for k in range(nvrt))
            f.write(" ".join(row) + " " + zvals + "\n")


# ----------------------------
# Plotting (MATLAB-like)
# ----------------------------
def plot_master(ax, m, hsm, zcor_m_masked):
    ax.plot(m, zcor_m_masked, "k-", linewidth=0.8)
    ax.plot(m, -hsm, "r.", markersize=6)
    nvrt_m = zcor_m_masked.shape[1]
    for i in range(len(m)):
        ax.plot(np.full(nvrt_m, m[i]), zcor_m_masked[i, :], "k", linewidth=0.6)
    ax.set_title("Master grid")
    ax.set_xlabel("Grid #")
    ax.set_ylabel("z (m)")
    ax.grid(True, alpha=0.25)


def plot_transect(ax, s, dp_line, zcor, title):
    ax.plot(s, zcor, "k-", linewidth=0.8)
    ax.plot(s, -dp_line, "r.", markersize=4)
    nvrt = zcor.shape[1]
    for i in range(len(s)):
        ax.plot(np.full(nvrt, s[i]), zcor[i, :], "k", linewidth=0.35)
    ax.set_title(title)
    ax.set_xlabel("Along transect distance (m)")
    ax.set_ylabel("z (m)")
    ax.grid(True, alpha=0.25)


def load_transect1(fname: str | Path = "transect1.out"):
    """
    Load transect1.out written in the same format as gen_vqs.f90:

      i, kbp, x, y, transect_len, dp, znd(:,nd)

    Column layout:
      col 0 : i
      col 1 : kbp
      col 2 : x
      col 3 : y
      col 4 : s (along-transect distance)
      col 5 : dp
      col 6+: z-coordinates (length = nvrt)

    Returns:
      kbp : (npbp,) int
      s   : (npbp,) float
      dp  : (npbp,) float
      zcor: (npbp, nvrt) float
    """
    arr = np.loadtxt(fname)

    kbp = arr[:, 1].astype(int)
    s   = arr[:, 4].astype(float)
    dp  = arr[:, 5].astype(float)
    zcor = arr[:, 6:].astype(float)

    return kbp, s, dp, zcor


def plot_VQS(
    vgrid_master: str = "vgrid_master.out",
    hgrid: str = "hgrid.gr3",
    vgrid: str = "vgrid.in",
    transect_bp: str = "transect.bp",
    transect1_out: str = "transect1.out",
    save_png: str = "vqs_diagnostics.png",
    show: bool = True,
):
    # ---- Load master and plot it ----
    m, nv_m, hsm, z_mas, zcor_m_masked = load_vgrid_master(vgrid_master)

    fig = plt.figure(figsize=(11, 8))
    ax1 = fig.add_subplot(2, 1, 1)
    plot_master(ax1, m, hsm, zcor_m_masked)

    ax2 = fig.add_subplot(2, 1, 2)

    bp_path = Path(transect_bp)
    if bp_path.exists():
        # ---- Build transect1.out from transect.bp + hgrid + vgrid ----
        ne, np_, x, y, dp, _, _ = read_hgrid_gr3(hgrid)
        nvrt, sigma_node, kbp = read_vgrid_in_ivcor1(vgrid)
        znd = sigma_to_znd(nvrt, sigma_node, dp, eta=0.0)

        xybp = read_transect_bp(transect_bp)
        s = cumulative_distance(xybp)

        imap = nearest_nodes(x, y, xybp)
        write_transect1_out(transect1_out, xybp, s, imap, dp, kbp, znd)

        # ---- Load the newly written transect1.out and plot ----
        kbp_t, s2, dp2, zcor1 = load_transect1(transect1_out)
        plot_transect(ax2, s2, dp2, zcor1, "Transect from transect.bp (nearest-node)")

    else:
        ax2.text(
            0.5, 0.5,
            f"No {transect_bp} found.\nProvide it to generate {transect1_out}.",
            ha="center", va="center", transform=ax2.transAxes
        )
        ax2.set_axis_off()

    fig.tight_layout()
    fig.savefig(save_png, dpi=200)
    print(f"Saved: {save_png}")
    if bp_path.exists():
        print(f"Wrote: {transect1_out}")
    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    plot_VQS(
        vgrid_master="/sciclone/schism10/feiye/STOFS3D-v8/I30/Vgrid/vgrid_master.out",
        hgrid="/sciclone/schism10/feiye/STOFS3D-v8/I30/Vgrid/hgrid.gr3",
        vgrid="/sciclone/schism10/feiye/STOFS3D-v8/I30/Vgrid/vgrid.in.old",
        transect_bp="/sciclone/schism10/feiye/STOFS3D-v8/I30/Vgrid0/transect.bp",
        transect1_out="/sciclone/schism10/feiye/STOFS3D-v8/I30/Vgrid/transect1.out",
        save_png="/sciclone/schism10/feiye/STOFS3D-v8/I30/Vgrid/vqs_diagnostics.png",
        show=True,
    )
