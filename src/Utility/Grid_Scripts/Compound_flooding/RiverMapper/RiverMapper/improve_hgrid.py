# %%
from pylib import schism_grid, grd2sms, proj_pts, read_schism_bpfile
# pylib is a python library that handles many schism-related file manipulations by Dr. Zhengui Wang
# , which can be installed by "pip install git+https://github.com/wzhengui/pylibs.git"
import shutil
import os
import numpy as np
from RiverMapper import improve_hgrid as original_file
from RiverMapper.Hgrid_extended import find_nearest_nd, hgrid_basic, read_schism_hgrid_cached, get_inp, propogate_nd, compute_ie_area
from RiverMapper.SMS import SMS_MAP
import pathlib
import copy
import subprocess
from glob import glob


iDiagnosticOutputs = False

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

def spring(gd):
    return gd

def reduce_bad_elements(
    gd=None, fixed_eles_id=np.array([], dtype=int), fixed_points_id=np.array([], dtype=int),
    area_threshold=150, skewness_threshold=10, nmax=6, output_fname=None
):

    n = 0
    while True:
        n += 1
        if n > nmax:
            break

        reduce_type = 2 + n % 2  #alternating two methods: (3) reduce 3 points of a triangle to one point; (2): reduce 2 points of a triangles' side
        print(f'\n>>>>Element reduction Iteration {n}')

        # *******************************************
        # <<Preparation>>
        # *******************************************
        # make a copy of the grid in case the current iteration fails
        gd0 = copy.deepcopy(gd)

        # calculate geometries; needs to update for every iteration
        gd.compute_bnd()
        print('computing additional geometries')
        gd.compute_ic3()
        gd.compute_area()
        gd.compute_nne()
        gd.compute_side(fmt=2)


        # *******************************************
        # Set fixed points and elements, i.e., those that cannot be changed
        # *******************************************
        # append boundary points to fixed_points fixed points
        if n > 1:  # clear input fixed points after the first iteration
            fixed_points_id = np.array([], dtype=int)
        fixed_points_id = np.sort(np.unique(np.array([*fixed_points_id, *gd.bndinfo.ip])))

        i_fixed_eles = np.zeros((gd.ne, ), dtype=bool)
        if len(fixed_points_id) > 0:
            fixed_eles_id = np.unique(gd.ine[fixed_points_id].flatten())
            fixed_eles_id = fixed_eles_id[fixed_eles_id>=0]
            i_fixed_eles[fixed_eles_id] = True

        # set quad nodes as fixed points
        # this is necessary because quad nodes may be moved when processing neighboring triangles,
        # even when the quad element itself set as fixed.
        quad_nodes = np.unique(gd.elnode[gd.i34==4, :])
        quad_nodes = quad_nodes[quad_nodes>=0]
        fixed_points_id = np.sort(np.unique(np.r_[fixed_points_id, quad_nodes]))

        i_fixed_points = np.zeros((gd.np, ), dtype=bool)
        i_fixed_points[fixed_points_id] = True

        # *******************************************
        # query grid quality; higher standards than in 'quality_check_hgrid'
        # *******************************************
        is_target_triangle = (gd.area < area_threshold) * (gd.i34 == 3)

        small_elnode_x = gd.x[gd.elnode[is_target_triangle, :3]]
        small_elnode_y = gd.y[gd.elnode[is_target_triangle, :3]]
        small_elnode_xy = small_elnode_x[:, :3] + 1j*small_elnode_y[:, :3]

        small_elnode_flatten = gd.elnode[is_target_triangle, :3].flatten()
        _, idx, counts = np.unique(small_elnode_xy.flatten(), return_counts=True, return_index=True)
        degenerate_count = np.zeros((gd.np, ), dtype=int)
        degenerate_count[small_elnode_flatten[idx]] = counts

        distj_padded = np.r_[gd.distj, np.nan]
        element_dl = distj_padded[gd.elside]
        aspect_ratios =  -np.log(abs(0.5 - np.nanmax(element_dl, axis=1) / (np.nansum(element_dl, axis=1) + 1e-12)))

        if reduce_type == 3:
            is_target_triangle = is_target_triangle * (aspect_ratios < skewness_threshold) * (~i_fixed_eles)
        else:
            is_target_triangle = is_target_triangle * (aspect_ratios >= skewness_threshold) * (~i_fixed_eles)
        target_triangle_id = np.argwhere(is_target_triangle)

        # *******************************************
        # <<loop control based on current grid quality>>
        # *******************************************
        if not is_target_triangle.any():
            print(f'No target elements found. Done fixing bad elements.')
            break

        # *******************************************
        # <<fix bad elements for this iteration>>
        # *******************************************
        print(f'trying to tweak {len(target_triangle_id)} Type-{reduce_type} elements\n')
        print(f'min element area = {np.sort(gd.area)[:20]} \n')

        gd_xy = gd.x + 1j*gd.y

        i_ele_shorted = np.zeros((gd.ne, ), dtype=bool)  # shorted elements, i.e., those not processed in the current iteration

        idx = np.argsort(gd.area[target_triangle_id].flatten())
        target_triangle_id_sorted = target_triangle_id[idx].flatten()

        for i, ie in enumerate(target_triangle_id_sorted):
            # progress bar
            if (len(target_triangle_id_sorted) > 100):
                if (i % int(len(target_triangle_id_sorted)/100) == 0 or i+1 == len(target_triangle_id_sorted)):
                    printProgressBar(i + 1, len(target_triangle_id_sorted), prefix = 'Progress:', suffix = 'Complete', length = 50)

            if i_ele_shorted[ie]:  # some elements are shorted out to prevent simultanesouly tweaking adjacent elements, or to prevent changing fixed points
                continue

            iee = ie
            for j in range(2):
                iee = np.unique(gd.ine[gd.elnode[iee]].flatten())
                iee = iee[iee>=0]
            i_ele_shorted[iee] = True

            # find degenerate triangles or sides
            if reduce_type == 3:
                ic3 = np.unique(gd.ic3[ie, :]); ic3 = ic3[ic3>=0]
                if (gd.i34[ic3] == 4).any():  # quad
                    i_ele_shorted[ie] = True  # skip
                    continue
                this_sorted_nodes = np.sort(gd.elnode[ie, :3].flatten())
                this_degenerate_count = degenerate_count[this_sorted_nodes]
                i_degenerate_node = np.zeros((3, ), dtype=bool)
                iee = [x for x in iee if not (x in gd.ic3[ie] or x==ie)]
            elif reduce_type == 2:  # find the shortest side as the degenerate side
                elsides = gd.elside[ie]
                dl = gd.distj[elsides[elsides>=0]]
                degenerate_side = np.argmin(dl)
                this_sorted_nodes = np.sort(gd.isidenode[elsides[degenerate_side], :])
                this_degenerate_count = degenerate_count[this_sorted_nodes]
                i_degenerate_node = np.zeros((2, ), dtype=bool)
                if gd.i34[gd.ic3[ie, degenerate_side]] == 4:  # quad
                    i_ele_shorted[ie] = True  # skip
                    continue
                iee = [x for x in iee if not (x == gd.ic3[ie, degenerate_side] or x==ie)]

            i_degenerate_node[np.argmax(this_degenerate_count)] = True
            move_node_id = this_sorted_nodes[~i_degenerate_node]
            if i_fixed_points[move_node_id].any():  # trying to move fixed nodes
                i_ele_shorted[ie] = True  # skip
                continue

            moved_nodes_id = this_sorted_nodes[~i_degenerate_node]
            degenerate_node_id = this_sorted_nodes[i_degenerate_node]
            gd_xy0_local = copy.deepcopy(gd_xy[moved_nodes_id])
            gd_xy[moved_nodes_id] = gd.x[degenerate_node_id] + 1j*gd.y[degenerate_node_id]  # reduce the last two points to the first point
            gd.x[moved_nodes_id] = np.real(gd_xy[moved_nodes_id])
            gd.y[moved_nodes_id] = np.imag(gd_xy[moved_nodes_id])
            area = compute_ie_area(gd, iee)
            if min(area) < 1e-4:  # revert the current edit
                i_ele_shorted[iee] = False
                i_ele_shorted[ie] = True
                gd_xy[moved_nodes_id] = gd_xy0_local

        # end loop ie in target_triangle_id_sorted

        # *******************************************
        # <<update basic grid info for output>>
        # *******************************************
        print(f'updating grid info ...')
        # update new node sequence
        gd_xy_unique, inv = np.unique(gd_xy, return_inverse=True)
        # update nodes
        gd.x = np.real(gd_xy_unique)
        gd.y = np.imag(gd_xy_unique)
        gd.dp = np.zeros((len(gd.x),), dtype=float)

        # update elements
        inv = np.r_[inv, -2, -2]  # pad two elements because elnode uses -2 as null
        # map node ids in ele map to new node sequence
        gd.elnode = inv[gd.elnode]

        # identify degenerate trianlges with duplicated node ids (new sequence)
        duplicate_nd_eles = np.zeros((gd.elnode.shape[0], ), dtype=bool)
        duplicate_nd_eles += (gd.elnode[:, 0] == gd.elnode[:, 1])
        duplicate_nd_eles += (gd.elnode[:, 0] == gd.elnode[:, 2])
        duplicate_nd_eles += (gd.elnode[:, 0] == gd.elnode[:, 3])
        duplicate_nd_eles += (gd.elnode[:, 1] == gd.elnode[:, 2])
        duplicate_nd_eles += (gd.elnode[:, 1] == gd.elnode[:, 3])
        duplicate_nd_eles += (gd.elnode[:, 2] == gd.elnode[:, 3])

        # remove degenerate ele
        gd.elnode = gd.elnode[~duplicate_nd_eles]
        gd.i34 = gd.i34[~duplicate_nd_eles]

        gd.ne = gd.elnode.shape[0]
        gd.np = gd.x.shape[0]

        # *******************************************
        # Finalize
        # *******************************************
        gd.compute_area()
        if min(gd.area) < 0.0:
            print(f'found negative elements: {np.argwhere(gd.area<0.0)+1}')
            print(f'reverting to previous grid ...')
            if iDiagnosticOutputs:
                grd2sms(gd, f'{os.path.dirname(output_fname)}/{pathlib.Path(output_fname).stem}.fix{n}.failed.2dm')
                grd2sms(gd0, f'{os.path.dirname(output_fname)}/{pathlib.Path(output_fname).stem}.fix{n}.pre-failed.2dm')
                print(f'final element areas: {np.sort(gd.area)}')
            gd = copy.deepcopy(gd0)
        else:
            if iDiagnosticOutputs:
                grd2sms(gd, f'{os.path.dirname(output_fname)}/{pathlib.Path(output_fname).stem}.fix{n}.2dm')
            gd = hgrid_basic(gd)
            pass

    # end of while loop

    gd = hgrid_basic(gd)
    if output_fname is not None:
        if pathlib.Path(output_fname).suffix == '.2dm':
            grd2sms(gd, output_fname)
        elif pathlib.Path(output_fname).suffix in ['.gr3', 'll']:
            gd.save(output_fname)
        else:
            raise Exception('suffix of output grid file name not supported.')

    return gd

def grid_element_relax(gd, target_points=None, niter=3, ntier=0, max_dist=50, min_area_allowed=1e-3, wdir=None, output_fname=None):
    if not os.path.exists(wdir):
        raise Exception(f'wdir: {wdir} not found')
    else:
        output_files = glob(f'{wdir}/*spring*gr3')
        if len(output_files)>0:
            for file in output_files:
                os.remove(file)

    # set fixed points (which won't be moved)
    if target_points is not None:
        if target_points.dtype == bool:
            pass
        elif target_points.dtype == int:
            tmp = np.zeros((gd.np, ), dtype=bool)
            tmp[target_points] = True
            target_points = tmp
        else:
            raise Exception('unsupported fixe_points type')
    else:  # all points can be moved if target_points are not specified
        target_points = np.ones((gd.np, ), dtype=bool)

    # fix boundary points too
    interior_points = np.ones((gd.np, ), dtype=bool)
    gd.compute_bnd()
    interior_points[gd.bndinfo.ip] = False

    print('preparing inputs')

    # find inp
    if not hasattr(gd, 'ine'):
        gd.compute_nne()

    i_relax = target_points * interior_points
    relax_points = np.argwhere(i_relax).flatten()

    # # Lloyd
    #     field = lloyd.Field(np.c_[gd.x[nd_nei], gd.y[nd_nei]])
    #     field.relax()
    #     new_xy = field.get_points()

    # Optional: multiple tiers (limited by max_dist)
    inp2 = get_inp(gd, ntiers=ntier).reshape(gd.np, -1)
    gd_xy = np.r_[gd.x + 1j*gd.y, 1e10+1j*1e10]
    nd_ids = np.array(range(gd.np))
    inp_distance = abs(gd_xy[inp2] - np.repeat(gd_xy[nd_ids].reshape(-1, 1), inp2.shape[1], axis=1))
    inp_shortrange = copy.deepcopy(inp2)
    inp_shortrange[inp_distance>max_dist] = -1

    relax_points = np.unique(inp_shortrange[target_points].flatten())
    relax_points = relax_points[relax_points>=0]
    i_relax = np.zeros((gd.np, ), dtype=bool)
    i_relax[relax_points] = True
    i_relax = i_relax * interior_points
    relax_points = np.argwhere(i_relax)

    ifixed = (~i_relax).astype(int)
    gd.write_hgrid(f'{wdir}/fixed.gr3', value=ifixed, fmt=1)

    # springing
    script_dir =  os.path.dirname(original_file.__file__)
    print(f'running grid_spring with {niter} iteration(s)...')
    p = subprocess.Popen(f'{script_dir}/grid_spring', cwd=wdir, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    p.stdin.write(f'{niter}\n{min_area_allowed}\n'.encode()) #expects a bytes type object
    p.communicate()[0]
    p.stdin.close()

    print('writing output files')
    gd = schism_grid(f'{wdir}/hgrid.spring.gr3')
    if output_fname is not None:
        if pathlib.Path(output_fname).suffix == '.2dm':
            grd2sms(gd, output_fname)
        elif pathlib.Path(output_fname).suffix in ['.gr3', 'll']:
            shutil.move(f'{wdir}/hgrid.spring.gr3', output_fname)
        else:
            raise Exception('suffix of output grid file name not supported.')

    gd.compute_area()
    print(f'Sorted area after springing: {np.sort(gd.area)[:50]}')

    return gd

def quality_check_hgrid(gd, outdir='./', small_ele_limit=5.0, skew_ele_minangle=0.5):
    '''
    Check several types of grid issues that may crash the model:
    Input hgrid needs to have meter unit
    '''

    print('\n-----------------quality check>')

    gd.compute_area()
    gd.compute_ctr()
    gd.compute_nne()

    # small and skew elements
    sorted_idx = np.argsort(gd.area)
    # small elements
    i_small_ele = gd.area < small_ele_limit
    n_small_ele = sum(i_small_ele)
    print(f'\n{n_small_ele} small (< {small_ele_limit} m2) elements:')
    small_ele = np.argwhere(i_small_ele).flatten()
    if n_small_ele > 0:
        print(f'Element id:            area,             xctr,             yctr')
        for ie in sorted_idx[:min(20, n_small_ele)]:
            print(f'Element {ie+1}: {gd.area[ie]}, {gd.xctr[ie]}, {gd.yctr[ie]}')

    # skew elements
    skew_ele = gd.check_skew_elems(angle_min=skew_ele_minangle, fmt=1, fname=None)
    print(f'\n{len(skew_ele)} skew (min angle < {skew_ele_minangle})')
    if len(skew_ele) > 0:
        sorted_idx = np.argsort(gd.area[skew_ele])
        sorted_skew_ele = np.sort(gd.area[skew_ele])
        print(f'Element id:            area,             xctr,             yctr')
        for i in sorted_idx[:min(20, len(sorted_skew_ele))]:
            ie = skew_ele[i]
            print(f'Element {ie+1}: {gd.area[ie]}, {gd.xctr[ie]}, {gd.yctr[ie]}')

    invalid = np.unique(np.array([*skew_ele, *small_ele]))
    if len(invalid) > 0:
        invalid_elnode = gd.elnode[invalid].reshape(-1, 1)
        invalid_elnode = np.unique(invalid_elnode)
        invalid_elnode = invalid_elnode[invalid_elnode>=0]
        invalid_neighbors = np.unique(gd.ine[invalid_elnode].reshape(-1, 1))
        invalid_neighbors = invalid_neighbors[invalid_neighbors>=0]
        i_invalid_nodes = np.zeros((gd.np, ), dtype=bool)
        i_invalid_nodes[invalid_elnode] = True  # set invalid nodes to be "not fixed", i.e., can be tweaked.

        # SMS_MAP(detached_nodes=np.c_[gd.xctr[invalid_neighbors], gd.yctr[invalid_neighbors], gd.yctr[invalid_neighbors]*0]).writer(f'{outdir}/invalid_element_relax.map')
    else:
        invalid_elnode = None
        invalid_neighbors = None
        i_invalid_nodes = None


    return {
        'hgrid': gd, 'invalid_nodes': invalid_elnode, 'invalid_elements': invalid_neighbors, 'i_invalid_nodes': i_invalid_nodes,
        'small_ele': small_ele, 'skew_ele': skew_ele
    }

def improve_hgrid(hgrid_name='', prj='esri:102008', load_bathy=False, nmax=3):
    '''
    Fix small and skew elements and bad quads
    prj: needs to specify hgrid's projection (the unit must be in meters)
    nmax: maximum number of rounds of fixing, most fixable elements can be fixed with nmax=4,
          nmax>4 ususally doesn't make additional improvements
    '''
    prj_name = prj.replace(':', '_')

    gd, dir_info = read_schism_hgrid_cached(hgrid_name, return_source_dir=True)
    dirname = dir_info['dir']
    file_basename = dir_info['basename']
    file_extension = dir_info['extension']

    # Fix invalid elements
    n_fix = 0
    while True:
        n_fix += 1

        grid_quality = quality_check_hgrid(gd, outdir=dirname)
        i_target_nodes = grid_quality['i_invalid_nodes']  # may be appended

        if i_target_nodes is None:  # all targets fixed
            print('\n -------------------------Done fixing invalid elements --------------------------------------------')
            break
        elif n_fix > nmax:  # maximum iteration reached, exit with leftovers
            print(' --------------------------------------------Done fixing invalid elements,')
            if len(grid_quality['small_ele']) > 0:
                print(f"Remaining small elements: {grid_quality['small_ele']}")
                SMS_MAP(detached_nodes=np.c_[gd.xctr[grid_quality['small_ele']], gd.yctr[grid_quality['small_ele']], gd.yctr[grid_quality['small_ele']]*0]).writer(f'{dirname}/small_ele.map')
            if len(grid_quality['skew_ele']) > 0:
                SMS_MAP(detached_nodes=np.c_[gd.xctr[grid_quality['small_ele']], gd.yctr[grid_quality['small_ele']], gd.yctr[grid_quality['small_ele']]*0]).writer(f'{dirname}/skew_ele.map')
                print(f"Remaining skew elements: {grid_quality['skew_ele']}")
            break
        else:  # fix targets

            print(f'\n----------------Fixing invalid elements, Round {n_fix}--------------------')

            # split bad quads
            print('\n ------------------- Splitting bad quads >')
            bp_name = f'{dirname}/bad_quad.bp'
            gd.check_quads(angle_min=60,angle_max=120,fname=bp_name)
            gd.split_quads(angle_min=60, angle_max=120)
            bad_quad_bp = read_schism_bpfile(fname=bp_name)
            print(f'{bad_quad_bp.nsta} bad quads split')
            if bad_quad_bp.nsta > 0:
                if iDiagnosticOutputs:
                    new_gr3_name = f"{dirname}/hgrid_split_quads.gr3"
                    gd.save(new_gr3_name)
                    print(f'the updated hgrid is saved as {new_gr3_name}')

                # quality check again since gd is updated
                i_target_nodes = quality_check_hgrid(gd, outdir=dirname)['i_invalid_nodes']
                if i_target_nodes is None: continue

            if n_fix == 1:
                # include intersection relax points only at Round 1
                inter_relax_pts = SMS_MAP(filename=f'{dirname}/intersection_relax.map').detached_nodes
                inter_relax_pts_x, inter_relax_pts_y = proj_pts(inter_relax_pts[:, 0], inter_relax_pts[:, 1], prj1='epsg:4326', prj2=prj)
                inter_relax_nd = find_nearest_nd(gd, np.c_[inter_relax_pts_x, inter_relax_pts_y])
                i_target_nodes[inter_relax_nd] = True

            print('\n ------------------- Reducing small/skew elements>')
            print(f'Number of target nodes: {sum(i_target_nodes)}')
            target_nodes_expand, i_target_nodes_expand = propogate_nd(gd, i_target_nodes, ntiers=3)
            gd = reduce_bad_elements(
                gd=gd, fixed_points_id=np.argwhere(~i_target_nodes_expand).flatten(),
                area_threshold=35,
                output_fname=f'{dirname}/{file_basename}_fix_bad_eles_round_{n_fix}.2dm'
            )

            # quality check again since gd is updated
            i_target_nodes = quality_check_hgrid(gd, outdir=dirname)['i_invalid_nodes']
            if i_target_nodes is None: continue

            if n_fix == 1:
                # re-find intersection nodes since gd is updated
                inter_relax_nd = find_nearest_nd(gd, np.c_[inter_relax_pts_x, inter_relax_pts_y])
                i_target_nodes[inter_relax_nd] = True


            print('\n ------------------- Relaxing remaining small/skew elements>')
            gd = grid_element_relax(
                gd=gd, target_points=i_target_nodes, niter=3, ntier=1, max_dist=20, wdir=dirname,
                output_fname=f'{dirname}/{file_basename}_relax_round_{n_fix}.2dm'
            )

            print(f'\n****************Done fixing invalid elements, Round {n_fix}*********************\n')

    grd2sms(gd, f'{dirname}/hgrid.2dm')
    gd.proj(prj0=prj, prj1='epsg:4326')
    gd.save(f'{dirname}/hgrid.ll')

    # load bathymetry
    # if load_bathy:
    #     os.chdir(dirname)
    #     os.system('python ./pload_depth.py')
    #     if not os.exist(f'{dirname}/hgrid.ll.new'):
    #         raise Exception('failed to load bathymetry')
    #     else:
    #         print('done loading bathymetry')

    pass

if __name__ == "__main__":
    # Sample usage
    gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I13r_LL/hgrid_xGEOID20b.gr3')
    gd.proj(prj0='epsg:4326', prj1='esri:102008')
    grd2sms(gd, '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ2/hgrid.102008.2dm')

    gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ2/hgrid.102008.fixed.2dm')
    quality_check_hgrid(gd)
    gd.save('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ2/hgrid.102008_fixed.gr3')

    gd = read_schism_hgrid_cached('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ2/hgrid.102008.fixed.gr3')
    gd.proj(prj0='esri:102008', prj1='epsg:4326')
    gd.save('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ2/hgrid.ll', fmt=1)

    # Sample usage
    improve_hgrid('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/V6_mesh_from_JZ2/v15.gr3')
