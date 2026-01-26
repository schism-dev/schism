#!/usr/bin/env python3
import numpy as np
from dataclasses import dataclass
from pathlib import Path

@dataclass
class HGrid:
    ne: int
    np: int
    x: np.ndarray      # (np,)
    y: np.ndarray      # (np,)
    dp: np.ndarray     # (np,) positive depth
    i34: np.ndarray    # (ne,)
    elnode: np.ndarray # (ne, 4) with 0 padding for triangles (or unused)

def read_hgrid_gr3(fname: str | Path) -> HGrid:
    fname = Path(fname)
    with fname.open("r") as f:
        _title = f.readline()
        ne, np_ = map(int, f.readline().split()[:2])

        x = np.empty(np_, dtype=float)
        y = np.empty(np_, dtype=float)
        dp = np.empty(np_, dtype=float)

        for i in range(np_):
            parts = f.readline().split()
            # j, x, y, dp
            x[i] = float(parts[1])
            y[i] = float(parts[2])
            dp[i] = float(parts[3])

        i34 = np.empty(ne, dtype=int)
        elnode = np.zeros((ne, 4), dtype=int)

        for e in range(ne):
            parts = f.readline().split()
            # j, i34, n1..n(i34)
            i34[e] = int(parts[1])
            nodes = list(map(int, parts[2:2 + i34[e]]))
            elnode[e, :i34[e]] = nodes  # 1-based node ids
    return HGrid(ne=ne, np=np_, x=x, y=y, dp=dp, i34=i34, elnode=elnode)


def build_master_vgrid2():
    """
    Alternative implementation of build_master_vgrid using more numpy features.
    """
    hsm_dict = {
        "1.0": 3, "2.0": 5, "3.0": 7,
        "4.0": 9, "6.0": 13, "8.0": 17, "12.0": 20, "18.0": 25,
        "25.0": 30, "33.0": 33, "42.0": 36, "52.0": 38, "67.0": 39,
        "83.0": 40, "100.0": 41, "150.0": 42, "230.0": 43, "350.0": 44,
        "1050.0": 53, "2000.0": 58, "5000.0": 63, "8500.0": 67,
    }
    m_vqs = len(hsm_dict)
    n_sd = 18
    dz_bot_min = 0.1
    a_vqs0 = 0.0
    etal = 0.0
    theta_b = 0.0

    hsm = np.array([float(k) for k in hsm_dict.keys()], dtype=float)
    nv_vqs = np.zeros(m_vqs, dtype=int)
    nv_vqs = np.array(list(hsm_dict.values()), dtype=int)

    if m_vqs < 2:
        raise ValueError("Check vgrid.in: m_vqs<2")
    if hsm[0] < 0:
        raise ValueError("hsm(1)<0")
    if np.any(np.diff(hsm) <= 0):
        raise ValueError("Check hsm: not strictly increasing")

    if etal <= -hsm[0]:
        raise ValueError("elev<hsm(1)")

    nvrt_m = int(nv_vqs[-1])
    z_mas = np.full((nvrt_m, m_vqs), -1.0e5, dtype=float)

    # Build first n_sd master grids using Dukhovskoy-like stretching
    for m in range(1, n_sd + 1):  # 1..n_sd in Fortran
        # theta_f logic (m is 1-based here)
        if m <= 7:
            theta_f = 0.0001
        elif m <= 17:
            theta_f = min(1.0, max(0.0001, (m - 4) / 10.0)) * 3.0
        else:
            theta_f = 4.4

        if m == 14: theta_f -= 0.1
        if m == 15: theta_f += 0.1
        if m == 16: theta_f += 0.55
        if m == 17: theta_f += 0.97

        nm = nv_vqs[m - 1]  # levels for this master
        for k in range(1, nm + 1):  # 1..nm
            sigma = (k - 1.0) / (1.0 - nm)  # same as Fortran
            # guard for extremely small theta_f
            sinh_tf = np.sinh(theta_f) if theta_f != 0 else 1.0
            cs = (1 - theta_b) * np.sinh(theta_f * sigma) / sinh_tf \
                 + theta_b * (np.tanh(theta_f * (sigma + 0.5)) - np.tanh(theta_f * 0.5)) / (2 * np.tanh(theta_f * 0.5))
            z_mas[k - 1, m - 1] = etal * (1 + sigma) + hsm[0] * sigma + (hsm[m - 1] - hsm[0]) * cs

    # Force downward steps for m=2..6 (Fortran)
    # Careful with indices: Fortran m=2..6 corresponds to 1..5 in 0-based
    for m in range(2, 7):  # 2..6 inclusive
        m0 = m - 1
        prev = m - 2
        nprev = nv_vqs[prev]
        nm = nv_vqs[m0]

        # z_mas(k,m)=min(z_mas(k,m),z_mas(k,m-1)) for k=1..nv_vqs(m-1)
        z_mas[:nprev, m0] = np.minimum(z_mas[:nprev, m0], z_mas[:nprev, prev])

        tmp = (z_mas[nm - 1, m0] - z_mas[nprev - 1, m0]) / (nm - nprev)
        for k in range(nprev + 1, nm + 1):  # k=nprev+1..nm
            z_mas[k - 1, m0] = z_mas[k - 2, m0] + tmp

    # Last 4 master grids copying + hard-coded deep z's
    # In Fortran: columns n_sd+1..n_sd+4 are 19..22 (1-based)
    col19 = n_sd      # 0-based index 18
    col20 = n_sd + 1  # 19
    col21 = n_sd + 2  # 20
    col22 = n_sd + 3  # 21

    # Copy top part from column n_sd (18, 1-based) => index n_sd-1 (17, 0-based)
    src = n_sd - 1
    z_mas[:nv_vqs[n_sd - 1], col19] = z_mas[:nv_vqs[n_sd - 1], src]
    z_mas[:nv_vqs[n_sd - 1], col20] = z_mas[:nv_vqs[n_sd - 1], src]
    z_mas[:nv_vqs[n_sd - 1], col21] = z_mas[:nv_vqs[n_sd - 1], src]
    z_mas[:nv_vqs[n_sd - 1], col22] = z_mas[:nv_vqs[n_sd - 1], src]

    # Fill deep parts (exact arrays from Fortran)
    def fill_deep(col, arr):
        start = nv_vqs[n_sd - 1]  # index where "-400 ..." begins (1+nv_vqs(n_sd) in Fortran)
        z_mas[start:start + len(arr), col] = np.array(arr, dtype=float)

    fill_deep(col19, [-400,-460,-520,-590,-660,-740,-830,-930,-1050])
    fill_deep(col20, [-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000])
    fill_deep(col21, [-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000,-2400,-2900,-3500,-4200,-5001])
    fill_deep(col22, [-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000,-2400,-2900,-3500,-4200,-5001,-6000,-7000,-8000,-9000])

    # Sanity: master nvrt
    nvrt_m = int(nv_vqs.max())

    return {
        "m_vqs": m_vqs,
        "n_sd": n_sd,
        "dz_bot_min": dz_bot_min,
        "a_vqs0": a_vqs0,
        "etal": etal,
        "hsm": hsm,
        "nv_vqs": nv_vqs,
        "nvrt_m": nvrt_m,
        "z_mas": z_mas,
    }


def build_master_vgrid():
    """
    Port of your master-grid section:
      m_vqs=22; n_sd=18; dz_bot_min=0.1; a_vqs0=0; etal=0
      hsm and nv_vqs arrays, plus z_mas generation and the "force downward steps" fix.

      n_sd is the number of standard masters (first n_sd columns); last 4 (m_vqs-n_sd) are custom deep grids.
    """
    m_vqs = 22
    n_sd = 18
    dz_bot_min = 0.1
    a_vqs0 = 0.0
    etal = 0.0
    theta_b = 0.0

    hsm = np.array([
        1., 2., 3.,
        4., 6., 8., 12., 18.,
        25., 33., 42., 52., 67.,
        83., 100., 150., 230., 350.,
        1050., 2000., 5000., 8500.
    ], dtype=float)
    assert len(hsm) == m_vqs

    nv_vqs = np.zeros(m_vqs, dtype=int)
    nv_vqs[:n_sd] = np.array([
        2, 3, 5,
        7, 9, 11, 13, 15,
        17, 17, 18, 18, 19,
        20, 21, 22, 24, 26
    ], dtype=int)
    nv_vqs[n_sd]     = nv_vqs[n_sd - 1] + 9   # n_sd+1 in Fortran
    nv_vqs[n_sd + 1] = nv_vqs[n_sd] + 5       # n_sd+2
    nv_vqs[n_sd + 2] = nv_vqs[n_sd + 1] + 5   # n_sd+3
    nv_vqs[n_sd + 3] = nv_vqs[n_sd + 2] + 4   # n_sd+4

    if m_vqs < 2:
        raise ValueError("Check vgrid.in: m_vqs<2")
    if hsm[0] < 0:
        raise ValueError("hsm(1)<0")
    if np.any(np.diff(hsm) <= 0):
        raise ValueError("Check hsm: not strictly increasing")

    if etal <= -hsm[0]:
        raise ValueError("elev<hsm(1)")

    nvrt_m = int(nv_vqs[-1])
    z_mas = np.full((nvrt_m, m_vqs), -1.0e5, dtype=float)

    # Build first n_sd master grids using Dukhovskoy-like stretching
    for m in range(1, n_sd + 1):  # 1..n_sd in Fortran
        # theta_f logic (m is 1-based here)
        if m <= 7:
            theta_f = 0.0001
        elif m <= 17:
            theta_f = min(1.0, max(0.0001, (m - 4) / 10.0)) * 3.0
        else:
            theta_f = 4.4

        if m == 14: theta_f -= 0.1
        if m == 15: theta_f += 0.1
        if m == 16: theta_f += 0.55
        if m == 17: theta_f += 0.97

        nm = nv_vqs[m - 1]  # levels for this master
        for k in range(1, nm + 1):  # 1..nm
            sigma = (k - 1.0) / (1.0 - nm)  # same as Fortran
            # guard for extremely small theta_f
            sinh_tf = np.sinh(theta_f) if theta_f != 0 else 1.0
            cs = (1 - theta_b) * np.sinh(theta_f * sigma) / sinh_tf \
                 + theta_b * (np.tanh(theta_f * (sigma + 0.5)) - np.tanh(theta_f * 0.5)) / (2 * np.tanh(theta_f * 0.5))
            z_mas[k - 1, m - 1] = etal * (1 + sigma) + hsm[0] * sigma + (hsm[m - 1] - hsm[0]) * cs

    # Force downward steps for m=2..6 (Fortran)
    # Careful with indices: Fortran m=2..6 corresponds to 1..5 in 0-based
    for m in range(2, 7):  # 2..6 inclusive
        m0 = m - 1
        prev = m - 2
        nprev = nv_vqs[prev]
        nm = nv_vqs[m0]

        # z_mas(k,m)=min(z_mas(k,m),z_mas(k,m-1)) for k=1..nv_vqs(m-1)
        z_mas[:nprev, m0] = np.minimum(z_mas[:nprev, m0], z_mas[:nprev, prev])

        tmp = (z_mas[nm - 1, m0] - z_mas[nprev - 1, m0]) / (nm - nprev)
        for k in range(nprev + 1, nm + 1):  # k=nprev+1..nm
            z_mas[k - 1, m0] = z_mas[k - 2, m0] + tmp

    # Last 4 master grids copying + hard-coded deep z's
    # In Fortran: columns n_sd+1..n_sd+4 are 19..22 (1-based)
    col19 = n_sd      # 0-based index 18
    col20 = n_sd + 1  # 19
    col21 = n_sd + 2  # 20
    col22 = n_sd + 3  # 21

    # Copy top part from column n_sd (18, 1-based) => index n_sd-1 (17, 0-based)
    src = n_sd - 1
    z_mas[:nv_vqs[n_sd - 1], col19] = z_mas[:nv_vqs[n_sd - 1], src]
    z_mas[:nv_vqs[n_sd - 1], col20] = z_mas[:nv_vqs[n_sd - 1], src]
    z_mas[:nv_vqs[n_sd - 1], col21] = z_mas[:nv_vqs[n_sd - 1], src]
    z_mas[:nv_vqs[n_sd - 1], col22] = z_mas[:nv_vqs[n_sd - 1], src]

    # Fill deep parts (exact arrays from Fortran)
    def fill_deep(col, arr):
        start = nv_vqs[n_sd - 1]  # index where "-400 ..." begins (1+nv_vqs(n_sd) in Fortran)
        z_mas[start:start + len(arr), col] = np.array(arr, dtype=float)

    fill_deep(col19, [-400,-460,-520,-590,-660,-740,-830,-930,-1050])
    fill_deep(col20, [-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000])
    fill_deep(col21, [-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000,-2400,-2900,-3500,-4200,-5001])
    fill_deep(col22, [-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000,-2400,-2900,-3500,-4200,-5001,-6000,-7000,-8000,-9000])

    # Sanity: master nvrt
    nvrt_m = int(nv_vqs.max())

    return {
        "m_vqs": m_vqs,
        "n_sd": n_sd,
        "dz_bot_min": dz_bot_min,
        "a_vqs0": a_vqs0,
        "etal": etal,
        "hsm": hsm,
        "nv_vqs": nv_vqs,
        "nvrt_m": nvrt_m,
        "z_mas": z_mas,
    }

def write_vgrid_master_out(fname: str | Path, hsm, nv_vqs, z_mas):
    fname = Path(fname)
    with fname.open("w") as f:
        for m in range(len(hsm)):
            nm = nv_vqs[m]
            # emulate Fortran: m, nv_vqs(m), hsm(m), z_mas(:,m)
            # only write up to nm levels (rest are -1e5)
            vals = " ".join(f"{z_mas[k, m]:12.4f}" for k in range(z_mas.shape[0]))
            f.write(f"{m+1:5d} {nm:5d} {hsm[m]:12.4f} {vals}\n")

def compute_vertical(hg: HGrid, master):
    """
    Port of:
      - dpmax check
      - per-node kbp and znd
      - sigma_vqs for shallow nodes computed immediately
      - later sigma_vqs computed from znd for deep nodes
    """
    hsm = master["hsm"]
    nv_vqs = master["nv_vqs"]
    z_mas = master["z_mas"]
    m_vqs = master["m_vqs"]
    dz_bot_min = master["dz_bot_min"]
    a_vqs0 = master["a_vqs0"]
    etal = master["etal"]

    dp = hg.dp
    np_ = hg.np

    dpmax = float(dp.max())
    if dpmax > float(hsm[m_vqs - 1]):
        raise ValueError(f"Max depth exceeds master depth: {dpmax} > {hsm[m_vqs-1]}")

    nvrt_m = master["nvrt_m"]
    eta2 = np.full(np_, etal, dtype=float)

    # Arrays (use 0-based internally; keep semantics)
    kbp = np.zeros(np_, dtype=int)                 # number of wet levels (like Fortran kbp(i))
    m0  = np.zeros(np_, dtype=int)                 # chosen master grid index (1-based in Fortran; we'll store 0-based)
    znd = np.full((nvrt_m, np_), -1.0e6, dtype=float)
    sigma_vqs = np.zeros((nvrt_m, np_), dtype=float)

    for i in range(np_):
        if dp[i] <= hsm[0]:  # shallow
            kbp_i = int(nv_vqs[0])
            kbp[i] = kbp_i
            nm = kbp_i
            for k in range(1, nm + 1):
                sigma = (k - 1.0) / (1.0 - nm)
                sig_t = a_vqs0 * sigma * sigma + (1.0 + a_vqs0) * sigma
                sigma_vqs[k - 1, i] = sig_t
                znd[k - 1, i] = sig_t * (eta2[i] + dp[i]) + eta2[i]
            continue

        # deep: find master bracket m such that hsm(m-1) < dp <= hsm(m) (Fortran m=2..m_vqs)
        found = False
        zrat = None
        mm = None
        for m in range(2, m_vqs + 1):  # 1-based
            if (dp[i] > hsm[m - 2]) and (dp[i] <= hsm[m - 1]):
                mm = m  # 1-based
                zrat = (dp[i] - hsm[m - 2]) / (hsm[m - 1] - hsm[m - 2])
                found = True
                break
        if not found or zrat is None:
            raise RuntimeError(f"Failed to find a master vgrid for node {i+1}, dp={dp[i]}")

        m0[i] = mm - 1  # store 0-based column for z_mas

        # Build znd until we hit bottom condition
        kbp_i = 0
        nm = int(nv_vqs[mm - 1])
        nprev = int(nv_vqs[mm - 2])  # mm>=2 always here
        z3_last = None

        for k in range(1, nm + 1):  # 1..nm
            # z1 = z_mas(min(k,nv_vqs(m-1)), m-1)
            kp = min(k, nprev)
            z1 = z_mas[kp - 1, mm - 2]
            z2 = z_mas[k - 1, mm - 1]
            z3 = z1 + (z2 - z1) * zrat
            z3_last = z3
            if z3 >= -dp[i] + dz_bot_min:
                znd[k - 1, i] = z3
            else:
                kbp_i = k
                break

        if kbp_i == 0:
            raise RuntimeError(
                f"Failed to find a bottom for node {i+1}, dp={dp[i]}; last z3={z3_last}"
            )

        # bottom level
        znd[kbp_i - 1, i] = -dp[i]
        kbp[i] = kbp_i

        # Check monotonic: znd(k-1) > znd(k) for k=2..kbp
        col = znd[:kbp_i, i]
        if np.any(col[:-1] <= col[1:]):
            raise RuntimeError(f"Inverted z at node {i+1}, dp={dp[i]}: {col}")

    # Extend beyond bottom for plotting (same as Fortran)
    for i in range(np_):
        znd[kbp[i]:, i] = -dp[i]

    return kbp, znd, sigma_vqs, eta2

def write_vgrid_in(fname: str | Path, hg: HGrid, master, kbp, znd, sigma_vqs, eta2):
    """
    Port of final vgrid.in writing (SELFE convention in your comment):
      ivcor=1
      nvrt = max(kbp)
      for each node i:
         write i, nvrt+1-kbp(i), sigma_vqs(kbp:1:-1,i)  (reversed, bottom->surface in file)
    """
    hsm = master["hsm"]
    nvrt = int(kbp.max())

    fname = Path(fname)
    with fname.open("w") as f:
        f.write("1\n")          # ivcor
        f.write(f"{nvrt}\n")    # nvrt

        for i in range(hg.np):
            if hg.dp[i] <= hsm[0]:
                # already has sigma_vqs from shallow computation
                pass
            else:
                # deep: reconstruct sigma from znd
                kb = int(kbp[i])
                sigma_vqs[0, i] = 0.0
                sigma_vqs[kb - 1, i] = -1.0
                denom = (eta2[i] + hg.dp[i])
                for k in range(2, kb):  # 2..kb-1 (1-based) => index 1..kb-2
                    sigma_vqs[k - 1, i] = (znd[k - 1, i] - eta2[i]) / denom

            # Check order: sigma(k) < sigma(k-1) for k=2..kbp (since sigma is negative downward)
            kb = int(kbp[i])
            col = sigma_vqs[:kb, i]
            if np.any(col[1:] >= col[:-1]):
                raise RuntimeError(f"Inverted sigma at node {i+1}, dp={hg.dp[i]}: {col}")

            # write: i, nvrt+1-kbp(i), sigma_vqs(kbp:1:-1)
            # sigma_vqs(kbp:1:-1) is bottom->surface in Fortran indexing
            sig_out = col[::-1]  # reverse
            n_dry = nvrt + 1 - kb
            f.write(f"{i+1:10d} {n_dry:10d} " + " ".join(f"{v:14.6f}" for v in sig_out) + "\n")

def write_nlev_gr3(fname: str | Path, hg: HGrid, kbp):
    fname = Path(fname)
    with fname.open("w") as f:
        f.write("# of levels at each node\n")
        f.write(f"{hg.ne} {hg.np}\n")
        for i in range(hg.np):
            f.write(f"{i+1} {hg.x[i]} {hg.y[i]} {int(kbp[i])}\n")
        for e in range(hg.ne):
            nn = hg.i34[e]
            nodes = hg.elnode[e, :nn]
            f.write(f"{e+1} {nn} " + " ".join(str(int(n)) for n in nodes) + "\n")


def gen_vqs(hgrid_file="hgrid.gr3", output_dir=None):
    """
    Generate vgrid.in and nlev.gr3 files for the STOFS3D model based on the provided hgrid.gr3 file.
    The generated files will be saved in the specified output directory (or current directory if None).
    """
    if output_dir is None:
        output_dir = Path.cwd()

    master = build_master_vgrid2()

    # Optional: save vgrid_master.out for debugging/plotting (as in Fortran)
    write_vgrid_master_out(f"{output_dir}/vgrid_master.out", master["hsm"], master["nv_vqs"], master["z_mas"])

    hg = read_hgrid_gr3(hgrid_file)
    kbp, znd, sigma_vqs, eta2 = compute_vertical(hg, master)

    print(f"Final nvrt = {kbp.max()}")
    # prisms count (same logic)
    # kbpl = max(kbp(elnode(1:i34(i),i))) and sum over elements
    # note: elnode is 1-based node ids -> convert to 0-based indices
    nprism = 0
    for e in range(hg.ne):
        nn = hg.i34[e]
        nids = hg.elnode[e, :nn] - 1
        kbpl = int(kbp[nids].max())
        nprism += kbpl
    print(f"# of prisms = {nprism}")
    print(f"Average # of layers = {nprism / hg.ne:.6f}")

    write_vgrid_in("vgrid.in", hg, master, kbp, znd, sigma_vqs, eta2)
    write_nlev_gr3("nlev.gr3", hg, kbp)


if __name__ == "__main__":
    gen_vqs()
