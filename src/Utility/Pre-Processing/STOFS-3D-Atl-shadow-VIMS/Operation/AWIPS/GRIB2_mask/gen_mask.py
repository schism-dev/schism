from pylib import schism_grid  # from ZG's pylib: pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pylibs4schism==0.1.10
from scipy import spatial
import numpy as np


if __name__ == "__main__":
    xyz = np.loadtxt('conus.xyz')
    xyz[xyz[:, 0]>180, 0] -= 360
    grid = schism_grid('./hgrid_split.gr3')

    ie, ip, acor = grid.compute_acor(np.c_[xyz[:, 0], xyz[:, 1]])
    ip[:] += -1  # node index start from 1
    ip[ie==-1, :] *= -1

    # fmt = '%12d '*5 + '%20.10f '*3
    # np.savetxt('/sciclone/schism10/feiye/STOFS3D-v4/AWIPS/AWIPS_mask_from_Yuji/mask0.txt', np.c_[np.array(range(len(ie)))+1, ie, ip, acor], fmt=fmt)
    fmt = '%12d '*4 + '%20.10f '*3
    np.savetxt('conus_mask2.txt', np.c_[np.array(range(len(ie)))+1, ip, acor], fmt=fmt)

    pass
