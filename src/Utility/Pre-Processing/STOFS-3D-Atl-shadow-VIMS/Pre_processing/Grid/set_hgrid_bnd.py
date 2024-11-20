import numpy as np
try:
    from pylib_experimental.schism_file import cread_schism_hgrid as schism_read
except ImportError:
    from pylib import read_schism_hgrid as schism_read


def convert_boundary_dict(boundary_dict):
    '''
    This function converts the boundary dictionary into
    the format for the boundary list input of pylib's hgrid class.
    For example:
        [[-94.21610701, 29.75730376], [-89.00481564, 30.39287886]]
    becomes
        [-94.21610701,-89.00481564,29.75730376,30.39287886 ]
    '''

    boundary_list = []

    for key in boundary_dict:
        boundary_points = np.array(boundary_dict[key])
        if boundary_points.shape != (2, 2):
            raise ValueError(
                'Each boundary should have two points' 
                ' (start and end, going counterclosewise along the boundary).'
            )
        boundary_list.append(
            [boundary_points[0, 0], boundary_points[1, 0],  # longitudes
                boundary_points[0, 1], boundary_points[1, 1]]  # latitudes
        )
    return boundary_list


def make_stofsv8_boundary(hgrid_obj, output_dir='./', write_hgrid=True):
    '''
    This function creates a boundary file for SCHISM.

    The first input is an hgird object (from pylib's read_schism_hgrid)

    The output is a boundary file, i.e., the boundary info segment of an hgrid.gr3
    '''

    boundary_dict = {
        "Atlantic": [[-58.93699999, 7.84699995], [-52.98670614, 46.75925769]],
        "Gulf of St. Lawrence": [[-55.39829625, 51.59437779], [-55.69225314, 52.09044046]],
        "St. Lawrence River": [[-71.20166881, 46.8421543], [-71.19458476, 46.81710242]],
    }

    boundary_list = convert_boundary_dict(boundary_dict)

    # set the boundary
    hgrid_obj.compute_bnd(boundary_list)
    hgrid_obj.write_bnd(f'{output_dir}/grd.bnd')
    if write_hgrid:
        print('Writing hgrid_with_bnd.gr3')
        hgrid_obj.save(f'{output_dir}/hgrid_with_bnd.gr3', fmt=1)  # fmt=1 for bnd


if __name__ == '__main__':
    wdir = '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v23.1'
    hg = schism_read(f'{wdir}/hgrid.gr3')
    make_stofsv8_boundary(hg, output_dir=wdir)
