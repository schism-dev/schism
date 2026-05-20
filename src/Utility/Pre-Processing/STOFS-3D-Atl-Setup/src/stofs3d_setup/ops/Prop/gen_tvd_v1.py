'''
Generate tvd.prop file for the hybrid WENO/ELM transport schemes in SCHISM

Rule:
    1 : WENO/TVD
    0 : upwind/ELM

Region convention:
    regions[0]              : nearshore region
    regions[1:]             : estuary regions or forced-upwind regions
    basename starts upwind  : forced tvd = 0
'''

import os
import socket
import numpy as np
import pylib

from pylib import read_schism_bpfile, inside_polygon
from pylib import schism_grid as schism_read


def _get_elnode_i34(hg):
    '''
    Ensure elnode shape = (ne, 4) and 0-based indexing
    '''
    elnode = np.asarray(hg.elnode).copy()
    i34 = np.asarray(hg.i34).astype(int).copy()

    # transpose if needed
    if elnode.shape[0] != hg.ne and elnode.shape[1] == hg.ne:
        elnode = elnode.T

    # convert 1-based -> 0-based if needed
    valid = elnode[elnode >= 0]
    if valid.size > 0:
        if valid.min() == 1 and valid.max() == hg.np:
            elnode = elnode - 1

    return elnode.astype(int), i34


def _build_node_to_elem(elnode, i34, np_node):
    '''
    Node-to-element connectivity
    '''
    node_to_elem = [[] for _ in range(np_node)]

    for ie in range(len(i34)):
        for nd in elnode[ie, :i34[ie]]:
            node_to_elem[nd].append(ie)

    return node_to_elem


def _build_side_neighbors(elnode, i34):
    '''
    Element side neighbors based on shared edges
    '''
    ne = len(i34)

    neighbors = [[] for _ in range(ne)]
    edge_owner = {}

    for ie in range(ne):
        nodes = list(elnode[ie, :i34[ie]])

        for k in range(i34[ie]):
            n1 = nodes[k]
            n2 = nodes[(k + 1) % i34[ie]]

            edge = tuple(sorted((n1, n2)))

            if edge in edge_owner:
                je = edge_owner[edge]

                neighbors[ie].append(je)
                neighbors[je].append(ie)

            else:
                edge_owner[edge] = ie

    return neighbors


def _split_regions(regions):
    '''
    Split input regions into:
        nearshore region
        estuary regions
        forced-upwind regions

    Any file whose basename starts with "upwind" is treated as forced-upwind.
    '''
    if len(regions) < 1:
        raise ValueError('regions must contain at least one nearshore region.')

    nearshore_region = regions[0]

    estuary_regions = []
    upwind_regions = []

    for region in regions[1:]:
        rname = os.path.basename(region)

        if rname.startswith('upwind'):
            upwind_regions.append(region)
        else:
            estuary_regions.append(region)

    return nearshore_region, estuary_regions, upwind_regions


def gen_tvd_prop(
    hg: pylib.schism_grid,
    regions: list,
    hmin_estu: float = 6.0,
    output_name: str = 'tvd.prop',
):
    '''
    Generate tvd.prop for hybrid WENO/ELM transport schemes

    Parameters
    ----------
    hg : schism_grid
        SCHISM grid object

    regions : list
        Region file list.

        regions[0] : nearshore region
        regions[1:] : estuary regions or forced-upwind regions

        Any region file whose basename starts with "upwind"
        is forced to tvd = 0.

    hmin_estu : float
        Cutoff depth for WENO/TVD

    output_name : str
        Output prop filename

    Returns
    -------
    None
    '''

    hg.compute_ctr()

    elnode, i34 = _get_elnode_i34(hg)

    node_to_elem = _build_node_to_elem(elnode, i34, hg.np)
    side_neighbors = _build_side_neighbors(elnode, i34)

    # -------------------------------------------------
    # Split regions
    # -------------------------------------------------
    nearshore_region, estuary_regions, upwind_regions = _split_regions(regions)

    print('Nearshore region:')
    print(f'  {nearshore_region}')

    print('Estuary regions:')
    for region in estuary_regions:
        print(f'  {region}')

    print('Forced upwind regions:')
    for region in upwind_regions:
        print(f'  {region}')

    # -------------------------------------------------
    # Nearshore region
    # -------------------------------------------------
    bp = read_schism_bpfile(nearshore_region, fmt=1)

    inear = inside_polygon(
        np.c_[hg.x, hg.y],
        bp.x,
        bp.y
    ).astype(int)

    # -------------------------------------------------
    # Estuary regions
    # -------------------------------------------------
    iest = np.zeros(hg.np, dtype=int)

    for region in estuary_regions:

        bp = read_schism_bpfile(region, fmt=1)

        idx = inside_polygon(
            np.c_[hg.x, hg.y],
            bp.x,
            bp.y
        ).astype(int)

        iest[idx == 1] = 1

    # -------------------------------------------------
    # Forced-upwind regions
    # -------------------------------------------------
    iupwind = np.zeros(hg.np, dtype=int)

    for region in upwind_regions:

        bp = read_schism_bpfile(region, fmt=1)

        idx = inside_polygon(
            np.c_[hg.x, hg.y],
            bp.x,
            bp.y
        ).astype(int)

        iupwind[idx == 1] = 1

    # -------------------------------------------------
    # Element flags
    # -------------------------------------------------
    inear_e = np.zeros(hg.ne, dtype=int)
    iest_e = np.zeros(hg.ne, dtype=int)
    iupwind_e = np.zeros(hg.ne, dtype=int)

    hmin_elem = np.zeros(hg.ne)

    for ie in range(hg.ne):

        nodes = elnode[ie, :i34[ie]]

        # if any node is inside nearshore region
        inear_e[ie] = np.max(inear[nodes])

        # only if all nodes are inside estuary region
        iest_e[ie] = np.min(iest[nodes])

        # if any node is inside forced-upwind region
        iupwind_e[ie] = np.max(iupwind[nodes])

        # minimum depth of the element
        hmin_elem[ie] = np.min(hg.dp[nodes])

    # -------------------------------------------------
    # Initial TVD/WENO flags
    #
    # 1 : WENO/TVD
    # 0 : upwind/ELM
    # -------------------------------------------------
    tvd = np.zeros(hg.ne, dtype=int)

    mask = (
        ((inear_e == 0) | (iest_e != 0))
        &
        (hmin_elem >= hmin_estu)
    )

    tvd[mask] = 1

    # -------------------------------------------------
    # Force upwind regions to tvd = 0
    # -------------------------------------------------
    tvd[iupwind_e == 1] = 0

    print(f'Initial WENO/TVD elements: {np.sum(tvd == 1)}')
    print(f'Initial upwind/ELM elements: {np.sum(tvd == 0)}')
    print(f'Forced-upwind elements: {np.sum(iupwind_e == 1)}')

    # -------------------------------------------------
    # Augment upwind zone by 1 extra layer offshore
    # -------------------------------------------------
    tvd0 = tvd.copy()

    for ie in range(hg.ne):

        if inear_e[ie] == 0 and tvd0[ie] == 0:

            nodes = elnode[ie, :i34[ie]]

            for nd in nodes:

                for je in node_to_elem[nd]:

                    tvd[je] = 0

    # Safety: forced-upwind regions remain tvd = 0
    tvd[iupwind_e == 1] = 0

    print(f'After 1-layer upwind augmentation:')
    print(f'  WENO/TVD elements: {np.sum(tvd == 1)}')
    print(f'  upwind/ELM elements: {np.sum(tvd == 0)}')

    # -------------------------------------------------
    # Remove isolated WENO elements
    # -------------------------------------------------
    while True:

        itouched = 0

        tvd0 = tvd.copy()

        for ie in range(hg.ne):

            if tvd0[ie] == 0:
                continue

            count_weno_neighbors = sum(
                tvd0[je] != 0
                for je in side_neighbors[ie]
            )

            if count_weno_neighbors <= 1:

                tvd[ie] = 0

                itouched += 1

        # Safety: forced-upwind regions remain tvd = 0
        tvd[iupwind_e == 1] = 0

        print(f'# of elem flipped = {itouched}')

        if itouched == 0:
            break

    # -------------------------------------------------
    # Final safety check
    # -------------------------------------------------
    tvd[iupwind_e == 1] = 0

    print('Final tvd.prop summary:')
    print(f'  WENO/TVD elements, tvd=1: {np.sum(tvd == 1)}')
    print(f'  upwind/ELM elements, tvd=0: {np.sum(tvd == 0)}')
    print(f'  forced-upwind elements: {np.sum(iupwind_e == 1)}')

    # -------------------------------------------------
    # Write output
    # -------------------------------------------------
    hg.write_prop(
        output_name,
        value=tvd,
        fmt='{:d}'
    )

    print(f'Wrote {output_name}')


def sample_usage():
    '''
    Sample usage
    '''

    gd = schism_read('hgrid.gr3')

    tvd_regions = [
        'iso_10m_edited.rgn',
        'Ches2.rgn',
        'DEBay.rgn',
        'Hudson.rgn',
        'upwind_Honduras.rgn',
        'upwind_east_Caribbean.rgn',
        'upwind_west_Caribbean.rgn',
    ]

    gen_tvd_prop(
        gd,
        regions=tvd_regions,
        hmin_estu=6.0,
        output_name='tvd.prop',
    )


if __name__ == '__main__':
    sample_usage()
