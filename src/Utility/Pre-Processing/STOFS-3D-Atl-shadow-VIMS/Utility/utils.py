import numpy as np
from matplotlib.tri import Triangulation

def split_quads(elnode=None):  # modified by FY
    '''
    Split quad elements to triangles and append additional elements to element table
    This script can be made much faster by using vector operation instead of the for-loop;
    just append additional elements to the end.
    '''
    if elnode is None:
        raise Exception('elnode should be a numpy array of (np,4)')

    elnode = np.ma.masked_values(elnode, -1) 
    if elnode.shape[1] == 4:
        eid = np.nonzero(~((np.isnan(elnode[:,-1]))|(elnode[:,-1]<0)))[0]
        elnode = np.r_[elnode[:,:3], np.c_[elnode[eid, 0][:, None], elnode[eid, 2:]]]
    return elnode.astype('int')

def triangulation(lon, lat, tris):
    if tris.max() >= len(lon): 
        tris -= 1
    return Triangulation(lon, lat, tris)
