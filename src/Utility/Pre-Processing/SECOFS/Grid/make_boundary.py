from pylib import read_schism_hgrid

def make_secofs_boundary(hgrid_obj, output_dir='./'):
    '''
    This function creates a boundary file for SCHISM.

    The first input is an hgird object (from pylib's read_schism_hgrid)
    
    The output is a boundary file, i.e., the boundary info segment of an hgrid.gr3
    '''

    boundary_dict = {
        "ocean_boundary": [
            [-86.76068604,-75.30909894,30.39105033,37.92420051],  # Atlantic Ocean
        ],
        "river_boundaries": [
            [-81.26937200,-81.26977200,32.52940200,32.52869300],  # Savannah River
            [-79.97328500,-79.97439500,33.21088100,33.21052800],  # Cooper River
        ]
    }

    # assemble the boundaries into a list
    boundary_list = []
    for key in boundary_dict:
        for boundary in boundary_dict[key]:
            boundary_list.append(boundary)

    # set the boundary
    hgrid_obj.compute_bnd(boundary_list)
    hgrid_obj.write_bnd(f'{output_dir}/grd.bnd')

if __name__ == '__main__':
    wdir = '/sciclone/home/hjyoo/schism10/task/task6_SECOFS/simulation/Whole_Domain/RUN06i_HJ/src/hgrid/'
    hgrid_obj = read_schism_hgrid(f'{wdir}/hgrid0.gr3')
    make_secofs_boundary(hgrid_obj, output_dir=wdir)