from pylib import schism_grid, sms2grd
import os
import pickle
import pathlib
import numpy as np
from scipy import spatial
from copy import deepcopy


def read_schism_hgrid_cached(gd_filename, overwrite_cache=False, return_source_dir=False):

    # gd_cache_fname = os.path.splitext(gd_filename)[0] + '.pkl'
    gd_cache_fname = gd_filename + '.pkl'

    if return_source_dir:
        dirname = os.path.dirname(gd_filename)
        file_basename = os.path.basename(gd_filename)
        file_extension = pathlib.Path(gd_filename).suffix


    if overwrite_cache or not os.path.exists(gd_cache_fname):
        if file_extension in ['.ll', '.gr3']:
            gd = schism_grid(gd_filename)
        elif file_extension == '.2dm':
            gd = sms2grd(gd_filename)

        with open(gd_cache_fname, 'wb') as file:
            pickle.dump(gd, file)
    else:
        with open(gd_cache_fname, 'rb') as file:
            gd = pickle.load(file)

    if return_source_dir:
        return gd, {"dir": dirname, "basename": file_basename, "extension": file_extension}
    else:
        return gd


def get_inp(gd, ntiers=1, return_grid=False):
    '''
    ntiers: number of tiers of neighboring nodes, ntiers=1 is the original inp,
            ntiers=0 is the original point
    '''
    if ntiers == 0:
        inp = np.array(range(gd.np))
    elif ntiers >=3:
        raise Exception('ntier>=3 is discouraged due to large memory use.')
    else:

        large_int = gd.np + 100
        
        if not hasattr(gd, 'ine'):
            gd.compute_nne()
        elnode_padded = np.r_[np.maximum(-1, gd.elnode), np.array([-1, -1, -1, -1]).reshape(1, 4)]
        inp = elnode_padded[gd.ine]
        inp = inp.reshape(inp.shape[0], inp.shape[1]*inp.shape[2])

        # squeeze
        inp[inp==-1] = large_int
        sorted = np.sort(inp, axis=1)
        diff = np.diff(sorted, axis=1)
        np.place(sorted[:, 1:], diff == 0, large_int)
        sorted.sort(axis=1)
        edge = np.argmax(sorted, axis=1).max()
        inp = sorted[:, :edge]
        inp[inp==large_int] = -1

        inp_padded = np.r_[inp, -np.ones((1, inp.shape[1]), dtype=int)]
        for i in range(ntiers-1):  # ntier=1 is the original inp, so no need for this loop
            inp = inp_padded[inp]
            inp = inp.reshape(inp.shape[0], inp.shape[1]*inp.shape[2])

        # futher processing
        inp[inp==-1] = large_int
        sorted = np.sort(inp, axis=1)
        diff = np.diff(sorted, axis=1)
        np.place(sorted[:, 1:], diff == 0, large_int)
        sorted.sort(axis=1)
        edge = np.argmax(sorted, axis=1).max()
        inp = sorted[:, :edge]
        inp[inp==large_int] = -1
    
    if return_grid:
        gd.inp = inp
        return inp, gd
    else:
        return inp

def propogate_nd(gd, nd_ids, ntiers=1):
    if ntiers < 1:
        raise Exception('no need to propogate')
    else:
        i_nd = np.zeros((gd.np,), dtype=bool)
        inp = get_inp(gd)
        for i in range(ntiers):
            nd_ids = np.unique(inp[nd_ids].flatten())
            nd_ids = nd_ids[nd_ids >= 0]
        i_nd[nd_ids] = True

    return nd_ids, i_nd
    
def get_bnd_nd_cached(gd, cache_file='./bnd_xy.pkl'):
    if not os.path.exists(cache_file):
        gd.compute_bnd()
        # cache boundary points coordinates
        bnd_x, bnd_y = gd.x[gd.bndinfo.ip], gd.y[gd.bndinfo.ip] 
        with open(cache_file, 'wb') as file:
            pickle.dump([bnd_x, bnd_y], file)
        bnd_nd = gd.bndinfo.ip
    else:
        # read boundary points coordinates from cache
        with open(cache_file, 'rb') as file:
            bnd_x, bnd_y = pickle.load(file)
        bnd_nd = spatial.cKDTree(np.c_[bnd_x, bnd_y]).query(gd.x, gd.y)[1]
    
    return bnd_nd

def hgrid_basic(gd):
    basic_attributes = ['source_file', 'ne', 'np', 'x', 'y', 'dp', 'i34', 'elnode', 'ns']
    gd_attributes = deepcopy(gd.__dict__)
    for attribute in gd_attributes:
        if attribute not in basic_attributes:
            delattr(gd, attribute)
    return gd

def compute_ie_area(gd, ie):
    fp=gd.elnode[ie,-1]<0;
    x1=gd.x[gd.elnode[ie,0]]; y1=gd.y[gd.elnode[ie,0]];
    x2=gd.x[gd.elnode[ie,1]]; y2=gd.y[gd.elnode[ie,1]];
    x3=gd.x[gd.elnode[ie,2]]; y3=gd.y[gd.elnode[ie,2]];
    x4=gd.x[gd.elnode[ie,3]]; y4=gd.y[gd.elnode[ie,3]]; x4[fp]=x1[fp]; y4[fp]=y1[fp]
    area=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)+(x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2
    return area

def find_nearest_nd(gd, pts):
    _, nd_ids = spatial.cKDTree(np.c_[gd.x, gd.y]).query(pts[:, :2])
    return nd_ids 


if __name__ == "__main__":
    gd = schism_grid('/sciclone/schism10/feiye/From_Nabi/Grid/02/hgrid.gr3')
    gd.compute_all()
    gd = hgrid_basic(gd)

    gd = sms2grd('/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/Relax_test4/Relax_test4.2dm')
    gd.split_quads()
    