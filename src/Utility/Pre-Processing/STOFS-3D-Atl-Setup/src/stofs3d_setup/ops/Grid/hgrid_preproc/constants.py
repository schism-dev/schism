"""Stable configuration used by the STOFS-3D hgrid workflow."""

STOFS_BOUNDARIES = {
    "Atlantic": [[-58.93699999, 7.84699995], [-52.98670614, 46.75925769]],
    "Gulf of St. Lawrence": [
        [-55.39829625, 51.59437779],
        [-55.69225314, 52.09044046],
    ],
    "St. Lawrence River": [
        [-71.20166881, 46.8421543],
        [-71.19458476, 46.81710242],
    ],
}

DUMMY_VGRID_LEVELS = 70
DUMMY_VGRID_MAX_DEPTH = 5000.0
DUMMY_VGRID_H_S = 5.0
DUMMY_VGRID_THETA_B = 0.75
DUMMY_VGRID_THETA_F = 4.0
