#!/usr/bin/env python
"""
This script fixes small and skew elements and bad quads in a SCHISM grid.

Usage:
Run in command line: python improve_hgrid.py <grid_file> --skewness_threshold <skewness_threshold> --area_threshold <area_threshold>
the grid file can be in 2dm, gr3, or ll format

Also see other sample usages in the main function.
"""

import socket
import sys
from pylib import schism_grid, sms2grd
from pylib import proj_pts, read_schism_bpfile, schism_bpfile

if 'gulf' in socket.gethostname():
    from pylib_experimental.schism_file import cread_schism_hgrid as read_schism_hgrid
    print('Using c++ function to accelerate hgrid reading')
else:
    from pylib import schism_grid as read_schism_hgrid
    print('Using python function to read hgrid')

# pylib is a python library that handles many schism-related file manipulations by Dr. Zhengui Wang
# , which can be installed by "pip install git+https://github.com/wzhengui/pylibs.git"
import shutil
import os
import numpy as np
from RiverMapper.Hgrid_extended import find_nearest_nd, hgrid_basic, get_inp, propogate_nd, compute_ie_area
from RiverMapper.SMS import SMS_MAP
import pathlib
import copy
import subprocess
from glob import glob
import argparse


iDiagnosticOutputs = False

def cmd_line_interface():
    parser = argparse.ArgumentParser(description="Input a grid file to fix small/skew elements.")

    parser.add_argument("grid", type=str, help="Path to the grid file to be fixed")
    parser.add_argument("--skewness_threshold", type=float, default=35.0, help="Maximum skewness allowed for an element, skewness assumes the definition in xmgredit, i.e., based on the ratio between the longest edge and the equivalent radius")
    parser.add_argument("--area_threshold", type=float, default=1.0, help="Minimum area allowed for an element")

    args = parser.parse_args()

    print("Grid file:", args.grid)
    print("Skewness_threshold:", args.skewness_threshold)
    print("Area threshold:", args.area_threshold)

    return args.grid, args.skewness_threshold, args.area_threshold

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

def reduce_bad_elements(
    gd=None, fixed_eles_id=np.array([], dtype=int), fixed_points_id=np.array([], dtype=int),
    area_threshold=None, skewness_threshold=None, nmax=4, output_fname=None
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
        print('computing geometries')
        gd.compute_all(fmt=1)
        gd.compute_bnd()

        # *******************************************
        # Set fixed points and elements, i.e., those that cannot be changed
        # *******************************************
        if n == 1:  # clear input fixed points after the first iteration
            fixed_points_id = np.sort(np.unique(np.array([*fixed_points_id])))
        else:
            fixed_points_id = np.array([], dtype=int)
        # append boundary points to fixed_points
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

        # More stringent than in quality_check_hgrid
        skew_ele = gd.check_skew_elems(threshold=skewness_threshold*0.8)

        print(f'trying to fix skew elements by edge swapping')
        i_swapped = np.zeros((gd.ne, ), dtype=bool)
        for i, ie in enumerate(skew_ele):

            if gd.i34[ie] == 4:  # quad
                continue

            isd = np.argmax(gd.distj[gd.elside[ie, :3]])
            nei = gd.ic3[ie, isd]

            # save the original element nodes and area
            area0, area_nei0 = compute_ie_area(gd, [ie, nei])
            elnode0 = copy.deepcopy(gd.elnode[ie, :3])
            elnode_nei0 = copy.deepcopy(gd.elnode[nei, :3])

            if i_swapped[ie] or i_swapped[nei]:  # already swapped
                continue

            if nei < 0 or gd.i34[nei] == 4:  # boundary or quad
                continue
            else:
                isd_nodes = gd.isidenode[gd.elside[ie, isd], :]  # nodes of the longest side
                non_share_node = elnode0[~np.isin(elnode0, isd_nodes)]
                if len(non_share_node) == 1:
                    non_share_node = non_share_node[0]
                else:
                    raise ValueError('non_share_node is not a scalar')
                non_share_node_nei = elnode_nei0[~np.isin(elnode_nei0, isd_nodes)]
                if len(non_share_node_nei) == 1:
                    non_share_node_nei = non_share_node_nei[0]
                else:
                    raise ValueError('non_share_node_nei is not a scalar')

                # re-arrange element nodes, starting from the non-share node
                idx = np.argwhere(gd.elnode[ie, :] == non_share_node).flatten()
                gd.elnode[ie, :3] = np.roll(gd.elnode[ie, :3], -idx)
                # swap the longest side
                gd.elnode[ie, 1] = non_share_node_nei

                # do the same for the neighbor
                idx = np.argwhere(gd.elnode[nei, :] == non_share_node_nei).flatten()
                gd.elnode[nei, :3] = np.roll(gd.elnode[nei, :3], -idx)
                gd.elnode[nei, 1] = non_share_node

            # check for improvement
            area, area_nei = compute_ie_area(gd, [ie, nei])
            if area < area0 or area_nei < area0:  # revert the swap
                gd.elnode[ie, :3] = elnode0
                gd.elnode[nei, :3] = elnode_nei0
            else:  # swap successful
                i_swapped[ie] = True
                i_swapped[nei] = True
                # print(f'swapped the common edge of elements {ie+1} and {nei+1}')

        print(f'swapped edges for {sum(i_swapped)} skew elements and their neighbors')

        # recalculate geometries
        gd.compute_all(fmt=1)

        # recalculate i_skew
        skew_ele = gd.check_skew_elems(threshold=skewness_threshold*0.8)
        i_skew = np.zeros((gd.ne, ), dtype=bool)
        i_skew[skew_ele] = True

        if reduce_type == 3:
            is_target_triangle = is_target_triangle * i_skew * (~i_fixed_eles)
        else:
            is_target_triangle = is_target_triangle * i_skew * (~i_fixed_eles)
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
                gd.grd2sms(f'{pathlib.Path(output_fname).parent}/{pathlib.Path(output_fname).stem}.fix{n}.failed.2dm')
                gd0.grd2sms(f'{pathlib.Path(output_fname).parent}/{pathlib.Path(output_fname).stem}.fix{n}.pre-failed.2dm')
                print(f'final element areas: {np.sort(gd.area)}')
            gd = copy.deepcopy(gd0)
        else:
            if iDiagnosticOutputs:
                gd.grd2sms(f'{pathlib.Path(output_fname).parent}/{pathlib.Path(output_fname).stem}.fix{n}.2dm')
            gd = hgrid_basic(gd)
            pass

    # end of while loop

    gd = hgrid_basic(gd)
    if output_fname is not None:
        if pathlib.Path(output_fname).suffix == '.2dm':
            gd.grd2sms(output_fname)
        elif pathlib.Path(output_fname).suffix in ['.gr3', 'll']:
            gd.save(output_fname)
        else:
            raise Exception('suffix of output grid file name not supported.')

    return gd


def grid_element_relax(gd, target_points=None, niter=3, ntier=0, max_dist=50, min_area_allowed=1e-3, wdir=None, output_fname=None):
    '''
    Relax specified element nodes to improve grid quality.
    The function prepares inputs for grid_spring, which is a Fortran executable,
    then calls grid_spring to do the actual work.

    Inputs:
        gd: schism_grid object
        target_points: a boolean array of size (gd.np, ) indicating whether a node is a target point
        niter: number of iterations
        ntier: number of tiers, i.e., number of layers of neighboring elements to be relaxed
        max_dist: maximum distance when searching for neighboring elements
        min_area_allowed: minimum area allowed for an element
        wdir: working directory of grid_spring (which is a Fortran executable)
        output_fname: output file name
    Outputs:
        gd: updated schism_grid object
        other diagnostic outputs are written to wdir
    '''

    if not os.path.exists(wdir):
        raise Exception(f'wdir: {wdir} not found')
    else:
        output_files = glob(f'{wdir}/*spring*gr3')
        if len(output_files)>0:
            for file in output_files:
                os.remove(file)

    # set target points, i.e., points that can be moved, as opposed to fixed points
    # user-specified target points
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

    # any interior points are potential target points, i.e., any boundary points are fixed
    interior_points = np.ones((gd.np, ), dtype=bool)
    gd.compute_bnd()
    interior_points[gd.bndinfo.ip] = False

    # relax points are user-specified target points that are also interior points
    i_relax = target_points * interior_points
    relax_points = np.argwhere(i_relax).flatten()

    print('preparing inputs')

    # find ine, i.e., neighboring elements of each node for all nodes
    if not hasattr(gd, 'ine'):
        gd.compute_nne()

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
    gd.write_hgrid(f'{wdir}/fixed.gr3', value=ifixed, fmt=0)

    # springing
    script_dir = os.path.dirname(__file__)
    print(f'running grid_spring with {niter} iteration(s)...')
    p = subprocess.Popen(f'{script_dir}/grid_spring', cwd=wdir, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    p.stdin.write(f'{niter}\n{min_area_allowed}\n'.encode())  # expects a bytes type object
    p.communicate()[0]
    p.stdin.close()

    print('writing output files')
    gd = schism_grid(f'{wdir}/hgrid.spring.gr3')
    if output_fname is not None:
        if pathlib.Path(output_fname).suffix == '.2dm':
            gd.grd2sms(output_fname)
        elif pathlib.Path(output_fname).suffix in ['.gr3', 'll']:
            shutil.move(f'{wdir}/hgrid.spring.gr3', output_fname)
        else:
            raise Exception('suffix of output grid file name not supported.')

    gd.compute_area()
    print(f'Sorted area after springing: {np.sort(gd.area)[:50]}')

    return gd


def cal_skewnewss(gd):
    '''
    Calculate skewness (longest_side_length/quivalent_element_radius) of each element
    '''
    if not hasattr(gd, 'dpe'):
        gd.compute_ctr()
    if not hasattr(gd, 'distj'):
        gd.compute_side(fmt=2)
    if not hasattr(gd, 'elside'):
        gd.compute_ic3()
    if not hasattr(gd, 'area'):
        gd.compute_area()

    distj = gd.distj[gd.elside]
    distj[gd.elside == -1] = 0  # side length

    gd.skewness = distj.max(axis=1)/np.sqrt(np.maximum(gd.area/np.pi, sys.float_info.epsilon))
    return gd.skewness


def quality_check_hgrid(gd, area_threshold=None, skewness_threshold=None):
    '''
    Check several types of grid issues that may crash the model,
    including small and skew elements, and negative area elements.
    The unit of the input hgrid needs to be in meters

    Inputs:
        gd: schism_grid object
        outdir: output directory
        area_threshold: minimum area allowed for an element
        skewness_threshold: maximum skewness allowed for an element,
            skewness assumes the definition in xmgredit, i.e.,
            based on the ratio between the longest edge and the equivalent radius
    '''

    print('\n-----------------quality check>')

    gd.compute_all(fmt=1)

    # negative area
    print(f'\n{sum(gd.area<0)} negative area elements.')
    # small and skew elements
    sorted_idx = np.argsort(gd.area)
    # small elements
    i_small_ele = gd.area < area_threshold
    n_small_ele = sum(i_small_ele)
    print(f'\n{n_small_ele} small (< {area_threshold} m2) elements:')
    small_ele = np.argwhere(i_small_ele).flatten()
    if n_small_ele > 0:
        print(f'Element id:            area,             xctr,             yctr')
        for ie in sorted_idx[:min(100, n_small_ele)]:
            print(f'Element {ie+1}: {gd.area[ie]}, {gd.xctr[ie]}, {gd.yctr[ie]}')

    # skew elements
    skew_ele = gd.check_skew_elems(threshold=skewness_threshold)
    print(f'\n{len(skew_ele)} skew (skewness >= {skewness_threshold})')
    if len(skew_ele) > 0:
        gd.skewness = cal_skewnewss(gd)
        sorted_idx = np.argsort(gd.skewness[skew_ele])[::-1]
        sorted_skew_ele = np.sort(gd.skewness[skew_ele])[::-1]
        print(f'Element id:            skewness,             xctr,             yctr')
        for i in sorted_idx[:min(100, len(sorted_skew_ele))]:
            ie = skew_ele[i]
            print(f'Element {ie+1}: {gd.skewness[ie]}, {gd.xctr[ie]}, {gd.yctr[ie]}')

    invalid = np.unique(np.array([*skew_ele, *small_ele]))
    if len(invalid) > 0:
        invalid_elnode = gd.elnode[invalid].reshape(-1, 1)
        invalid_elnode = np.unique(invalid_elnode)
        invalid_elnode = invalid_elnode[invalid_elnode>=0]
        invalid_neighbors = np.unique(gd.ine[invalid_elnode].reshape(-1, 1))
        invalid_neighbors = invalid_neighbors[invalid_neighbors>=0]
        i_invalid_nodes = np.zeros((gd.np, ), dtype=bool)
        i_invalid_nodes[invalid_elnode] = True  # set invalid nodes to be "not fixed", i.e., can be tweaked.
        # SMS_MAP(detached_nodes=np.c_[
        #     gd.xctr[invalid_neighbors], gd.yctr[invalid_neighbors], gd.yctr[invalid_neighbors]*0
        # ]).writer(f'{outdir}/invalid_element_relax.map')
    else:
        invalid_elnode = None
        invalid_neighbors = None
        i_invalid_nodes = None

    return {
        'hgrid': gd, 'invalid_nodes': invalid_elnode, 'invalid_elements': invalid_neighbors, 'i_invalid_nodes': i_invalid_nodes,
        'small_ele': small_ele, 'skew_ele': skew_ele
    }

def write_diagnostics(outdir=None, grid_quality=None, hgrid_ref=None):
    '''hgrid_ref: reference hgrid, used to write diagnostic outputs,
    , which can be in a different projection than the hgrid associated with grid_quality'''

    if not hasattr(hgrid_ref, 'xctr'):
        hgrid_ref.compute_ctr()

    print(f"Remaining small elements: {grid_quality['small_ele']}")
    if len(grid_quality['small_ele']) > 0:
        SMS_MAP(detached_nodes=np.c_[hgrid_ref.xctr[grid_quality['small_ele']], hgrid_ref.yctr[grid_quality['small_ele']], hgrid_ref.yctr[grid_quality['small_ele']]*0]).writer(f'{outdir}/small_ele.map')
        small_ele_bp = schism_bpfile(x=hgrid_ref.xctr[grid_quality['small_ele']], y=hgrid_ref.yctr[grid_quality['small_ele']])
        small_ele_bp.write(f'{outdir}/small_ele.bp')

    print(f"Remaining skew elements: {grid_quality['skew_ele']}")
    if len(grid_quality['skew_ele']) > 0:
        SMS_MAP(detached_nodes=np.c_[hgrid_ref.xctr[grid_quality['skew_ele']], hgrid_ref.yctr[grid_quality['skew_ele']], hgrid_ref.yctr[grid_quality['skew_ele']]*0]).writer(f'{outdir}/skew_ele.map')
        skew_ele_bp = schism_bpfile(x=hgrid_ref.xctr[grid_quality['skew_ele']], y=hgrid_ref.yctr[grid_quality['skew_ele']])
        skew_ele_bp.write(f'{outdir}/skew_ele.bp')


def improve_hgrid(gd, prj='esri:102008', skewness_threshold=35, area_threshold=1,n_intersection_fix=0, nmax=5):
    '''
    Fix small and skew elements and bad quads
    prj: needs to specify hgrid's projection (the unit must be in meters)
    nmax: maximum number of rounds of fixing, most fixable elements can be fixed with nmax=4,
          nmax>4 ususally doesn't make additional improvements
    n_intersection_fix: number of rounds of relaxing intersection points
          default is 0, i.e., not relaxing intersection points
          set a large number (e.g., 999) to relax intersection points in every round
    '''
    dirname = pathlib.Path(gd.source_file).resolve().parent
    file_basename = os.path.basename(gd.source_file)

    # Fix invalid elements
    n_fix = 0
    while True:
        n_fix += 1

        grid_quality = quality_check_hgrid(gd, area_threshold=area_threshold, skewness_threshold=skewness_threshold)
        i_target_nodes = grid_quality['i_invalid_nodes']  # may be appended

        if i_target_nodes is None:  # all targets fixed
            print('\n -------------------------Done fixing invalid elements --------------------------------------------')
            break
        elif n_fix > nmax:  # maximum iteration reached, exit with leftovers
            print(' ----------------------------- Done fixing invalid elements, but with leftovers ---------------------')
            write_diagnostics(outdir=dirname, grid_quality=grid_quality, hgrid_ref=gd)
            break
        else:  # fix targets

            print(f'\n----------------Fixing invalid elements, Round {n_fix}--------------------')

            if sum(gd.i34 == 4) > 0:  # only check quads if there are any
                print('\n ------------------- Splitting bad quads >')
                bp_name = f'{dirname}/bad_quad.bp'
                # split bad quads
                bad_quad_idx = gd.check_quads(angle_min=50,angle_max=130,fname=bp_name)
                if bad_quad_idx is not None and len(bad_quad_idx) > 0:
                    gd.split_quads(angle_min=50, angle_max=130)

                    # outputs from split quads
                    bad_quad_bp = read_schism_bpfile(fname=bp_name)
                    print(f'{bad_quad_bp.nsta} bad quads split')
                    if bad_quad_bp.nsta > 0:
                        if iDiagnosticOutputs:
                            new_gr3_name = f"{dirname}/hgrid_split_quads.gr3"
                            gd.save(new_gr3_name)
                            print(f'the updated hgrid is saved as {new_gr3_name}')
                        # quality check again since gd is updated
                        grid_quality = quality_check_hgrid(gd, area_threshold=area_threshold, skewness_threshold=skewness_threshold)['i_invalid_nodes']
                        if i_target_nodes is None: continue

            gd.compute_all(fmt=1)

            print('\n ----------------------------- include intersection relax points ----------------------')
            if n_fix <= n_intersection_fix:
                inter_relax_pts = SMS_MAP(filename=f'{dirname}/total_intersection_joints.map').detached_nodes
                inter_relax_pts_x, inter_relax_pts_y = proj_pts(inter_relax_pts[:, 0], inter_relax_pts[:, 1], prj1='epsg:4326', prj2=prj)
                inter_relax_nd = find_nearest_nd(gd, np.c_[inter_relax_pts_x, inter_relax_pts_y])
                i_target_nodes[inter_relax_nd] = True

            print('\n ------------------- Reducing small/skew elements>')
            print(f'Number of target nodes: {sum(i_target_nodes)}')
            target_nodes_expand, i_target_nodes_expand = propogate_nd(gd, i_target_nodes, ntiers=3)
            gd = reduce_bad_elements(
                gd=gd, fixed_points_id=np.argwhere(~i_target_nodes_expand).flatten(),
                skewness_threshold=skewness_threshold * 0.8,  # more stringent than in quality_check_hgrid
                area_threshold=area_threshold * 10,  # more stringent than in quality_check_hgrid
                output_fname=f'{dirname}/{file_basename}_fix_bad_eles_round_{n_fix}.2dm'
            )

            # quality check again since gd is updated
            i_target_nodes = quality_check_hgrid(gd, area_threshold=area_threshold, skewness_threshold=skewness_threshold)['i_invalid_nodes']
            if i_target_nodes is None: continue

            if n_fix <= n_intersection_fix:
                # re-find intersection nodes since gd is updated
                inter_relax_nd = find_nearest_nd(gd, np.c_[inter_relax_pts_x, inter_relax_pts_y])
                i_target_nodes[inter_relax_nd] = True


            print('\n ------------------- Relaxing remaining small/skew elements>')
            gd = grid_element_relax(
                gd=gd, target_points=i_target_nodes, niter=3, ntier=1, max_dist=20, wdir=dirname,
                output_fname=f'{dirname}/{file_basename}_relax_round_{n_fix}.2dm'
            )

            print(f'\n****************Done fixing invalid elements, Round {n_fix}*********************\n')
    # end while loop

    if prj != 'epsg:4326':
        gd_x, gd_y = gd.x, gd.y  # save a copy of the original coordinates
        gd.proj(prj0=prj, prj1='epsg:4326')
        gd_lon, gd_lat = gd.x, gd.y

    if sum(gd.i34 == 4) > 0:  # only check quads if there are any
        print('\n ------------------- Splitting bad quads in lon/lat >')
        bp_name = f'{dirname}/bad_quad.bp'
        # split bad quads
        gd.check_quads(angle_min=60,angle_max=120,fname=bp_name)
        gd.split_quads(angle_min=60, angle_max=120)
        # outputs from split quads
        bad_quad_bp = read_schism_bpfile(fname=bp_name)
        print(f'{bad_quad_bp.nsta} bad quads split')
        if bad_quad_bp.nsta > 0:
            if iDiagnosticOutputs:
                new_gr3_name = f"{dirname}/hgrid_split_quads.gr3"
                gd.save(new_gr3_name)
                print(f'the updated hgrid is saved as {new_gr3_name}')
            # quality check again since gd is updated
            if prj != 'epsg:4326':
                gd.x, gd.y = gd_x, gd_y  # revert to original coordinates
            i_target_nodes = quality_check_hgrid(gd, area_threshold=area_threshold, skewness_threshold=skewness_threshold)['i_invalid_nodes']

    print('\n ------------------- Outputting final hgrid >')

    return gd

    pass


def sample1():
    '''
    Sample usage 1 without command line interface.
    The input grid is in *.ll format and lon/lat
    '''

    grid_dir = '/sciclone/schism10/feiye/STOFS3D-v8/R15e_v7/'
    grid_file = f'{grid_dir}/hgrid_xy_transferred.ll'

    # this test may find any potential boundary issues
    gd = read_schism_hgrid(grid_file)
    gd.compute_area()
    gd.compute_bnd(method=1)

    # manually set parameters
    skewness_threshold = 44.5
    area_threshold = 5
    gd = schism_grid(grid_file)

    gd_ll = copy.deepcopy(gd)
    # gd_ll.proj(prj0='esri:102008', prj1='epsg:4326')  # reproject to lon/lat if necessary
    gd_meter = copy.deepcopy(gd)
    gd_meter.proj(prj0='epsg:4326', prj1='esri:102008')  # reproject to meters if necessary

    grid_quality = quality_check_hgrid(gd_meter, area_threshold=area_threshold, skewness_threshold=skewness_threshold)
    write_diagnostics(outdir=grid_dir, grid_quality=grid_quality, hgrid_ref=gd_ll)

    # improve grid quality
    gd = improve_hgrid(gd_meter, n_intersection_fix=0, area_threshold=area_threshold, skewness_threshold=skewness_threshold, nmax=4)

    gd.grd2sms(f'{grid_dir}/hgrid.2dm')
    gd.save(f'{grid_dir}/hgrid.ll', fmt=1)


def sample2():
    '''
    Sample usage 2 without command line interface
    The input grid is in *.2dm format and esri:102008 projection
    '''
    grid_dir = '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v54_s2v3/Improve/'
    grid_file = f'{grid_dir}/v54_s2v3.2dm'

    # read hgrid
    # gd = read_schism_hgrid(grid_file)  # esri:102008
    gd = sms2grd(grid_file)
    gd.source_file = grid_file
    gd_ll = copy.deepcopy(gd)  # save a copy

    # gd_ll.proj(prj0='esri:102008', prj1='epsg:4326')
    gd.proj(prj0='epsg:4326', prj1='esri:102008')  # reproject to meters if necessary

    # this test may find any potential boundary issues
    gd.compute_area()
    gd.compute_bnd(method=1)

    # manually set parameters
    skewness_threshold = 25
    area_threshold = 5

    grid_quality = quality_check_hgrid(gd, area_threshold=area_threshold, skewness_threshold=skewness_threshold)
    write_diagnostics(outdir=grid_dir, grid_quality=grid_quality, hgrid_ref=gd_ll)

    # improve grid quality
    improve_hgrid(gd, n_intersection_fix=0, area_threshold=area_threshold, skewness_threshold=skewness_threshold, nmax=2)

    # project back to lon/lat if necessary
    gd.proj(prj0='esri:102008', prj1='epsg:4326')
    gd.grd2sms(f'{grid_dir}/hgrid.2dm')
    gd.save(f'{grid_dir}/hgrid.ll', fmt=1)


def sample3():
    '''
    Sample usage 3 without command line interface.
    The input grid is in esri:102008
    '''

    grid_dir = '/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/v32d/Improve/'
    grid_file = f'{grid_dir}/v32d.2dm'

    # this test may find any potential boundary issues
    gd = sms2grd(grid_file)
    gd.source_file = grid_file
    gd.compute_area()
    gd.compute_bnd(method=1)

    # manually set parameters
    skewness_threshold = 100
    area_threshold = 1

    gd_meter = copy.deepcopy(gd)
    gd_ll = copy.deepcopy(gd)
    gd_ll.proj(prj0='esri:102008', prj1='epsg:4326')  # reproject to lon/lat if necessary

    grid_quality = quality_check_hgrid(gd_meter, area_threshold=area_threshold, skewness_threshold=skewness_threshold)
    write_diagnostics(outdir=grid_dir, grid_quality=grid_quality, hgrid_ref=gd_ll)

    # improve grid quality
    gd = improve_hgrid(gd_meter, n_intersection_fix=0, area_threshold=area_threshold, skewness_threshold=skewness_threshold, nmax=4)

    gd.grd2sms(f'{grid_dir}/hgrid.2dm')
    gd.save(f'{grid_dir}/hgrid.gr3', fmt=1)
    gd_ll = copy.deepcopy(gd)  # needed because gd is updated after improving grid quality
    gd_ll.proj(prj0='esri:102008', prj1='epsg:4326')  # reproject to lon/lat if necessary
    gd_ll.save(f'{grid_dir}/hgrid.ll', value=gd.dp, fmt=1)


def main():
    '''
    # Sample usage , the horizontal coordinate unit of hgrid must be in meters
    # if you grid is in lon/lat, convert it to meters first, e.g.:
    # gd.proj(prj0='epsg:4326', prj1='esri:102008')
    '''
    grid_file, skewness_threshold, area_threshold = cmd_line_interface()
    gd = read_schism_hgrid(grid_file)

    # sanity check for illegal boundaries, in case of which this step will hang
    gd.compute_bnd(method=1)

    # improve grid quality
    improve_hgrid(gd, n_intersection_fix=0, area_threshold=area_threshold, skewness_threshold=skewness_threshold, nmax=4)

    grid_dir = pathlib.Path(grid_file).resolve().parent
    if pathlib.Path(grid_file).suffix in ['.2dm']:
        gd.grd2sms(f'{grid_dir}/improved_{pathlib.Path(grid_file).name}')
    elif pathlib.Path(grid_file).suffix in ['.gr3', '.ll']:
        gd.save(f'{grid_dir}/improved_{pathlib.Path(grid_file).name}', fmt=1)


if __name__ == "__main__":
    sample3()
    # main()

