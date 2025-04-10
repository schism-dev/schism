'''
Usage: python generate_adcirc.py -h

For example:

(on WCOSS2) python generate_adcirc.py --input_filename t12z.fields.out2d_nowcast.nc --input_city_identifier_file ./Shapefiles/city_poly.node_id.txt --output_dir .
(on other clusters) python generate_adcirc.py --input_filename ./outputs/out2d_1.nc --input_city_identifier_file ./Shapefiles/city_poly.shp --output_dir ./extract/  --datum "xGEOID20B"

Outputs:
./extract/schout_adcirc_*.nc, which are in ADCIRC's format

'''
from time import time
import argparse
from pathlib import Path
import copy
import os
import errno

import matplotlib as mpl
import shapefile

import numpy as np
from netCDF4 import Dataset

def signa(x,y):
    '''
        compute signed area for triangles along the last dimension (x[...,0:3],y[...,0:3])
    '''
    if x.ndim==1:
        area=((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]))/2
    elif x.ndim==2:
        # area=((x[:,0]-x[:,2])*(y[:,1]-y[:,2])-(x[:,1]-x[:,2])*(y[:,0]-y[:,2]))/2
        area=((x[...,0]-x[...,2])*(y[...,1]-y[...,2])-(x[...,1]-x[...,2])*(y[...,0]-y[...,2]))/2
    area=np.squeeze(area)
    return area

def inside_polygon(pts,px,py,fmt=0,method=0):
    '''
    check whether points are inside polygons

    usage: sind=inside_polygon(pts,px,py):
       pts[npt,2]: xy of points
       px[npt] or px[npt,nploy]: x coordiations of polygons
       py[npt] or py[npt,nploy]: y coordiations of polygons
       (npt is number of points, nploy is number of polygons)

       fmt=0: return flags "index[npt,nploy]" for whether points are inside polygons (1 means Yes, 0 means No)
       fmt=1: only return the indices "index[npt]" of polygons that pts resides in
             (if a point is inside multiple polygons, only one indice is returned; -1 mean pt outside of all Polygons)

       method=0: use mpl.path.Path
       method=1: use ray method explicitly

    note: For method=1, the computation time is proportional to npt**2 of the polygons. If the geometry
          of polygons are too complex, dividing them to subregions will increase efficiency.
    '''
    #----use ray method-----------------------------
    #check dimension
    if px.ndim==1:
       px=px[:,None]; py=py[:,None]

    #get dimensions
    npt=pts.shape[0]; nv,npy=px.shape

    if nv==3 and fmt==1:
       #for triangles, and only return indices of polygons that points resides in
       px1=px.min(axis=0); px2=px.max(axis=0); py1=py.min(axis=0); py2=py.max(axis=0)

       sind=[];
       for i in np.arange(npt):
           pxi=pts[i,0]; pyi=pts[i,1]
           sindp=np.nonzero((pxi>=px1)*(pxi<=px2)*(pyi>=py1)*(pyi<=py2))[0]; npy=len(sindp)
           if npy==0:
               sind.append(-1)
           else:
               isum=np.ones(npy)
               for m in np.arange(nv):
                   xi=np.c_[np.ones(npy)*pxi,px[m,sindp],px[np.mod(m+1,nv),sindp]]
                   yi=np.c_[np.ones(npy)*pyi,py[m,sindp],py[np.mod(m+1,nv),sindp]]
                   area=np.signa(xi,yi)
                   fp=area<0; isum[fp]=0;
               sindi=np.nonzero(isum!=0)[0]

               if len(sindi)==0:
                   sind.append(-1)
               else:
                   sind.append(sindp[sindi[0]])
       sind=np.array(sind)

    else:
        if method==0:
            sind=[]
            for m in np.arange(npy):
                sindi=mpl.path.Path(np.c_[px[:,m],py[:,m]]).contains_points(pts)
                sind.append(sindi)
            sind=np.array(sind).T+0  #convert logical to int
        elif method==1  :
            #using ray method explicitly
            sind=np.ones([npt,npy])
            x1=pts[:,0][:,None]; y1=pts[:,1][:,None]
            for m in np.arange(nv):
                x2=px[m,:][None,:]; y2=py[m,:][None,:]; isum=np.zeros([npt,npy])
                # sign_x1_x2=sign(x1-x2)
                for n in np.arange(1,nv-1):
                    x3=px[(n+m)%nv,:][None,:]; y3=py[(n+m)%nv,:][None,:]
                    x4=px[(n+m+1)%nv,:][None,:]; y4=py[(n+m+1)%nv,:][None,:]

                    #formulation for a ray to intersect with a line
                    fp1=((y1-y3)*(x2-x1)+(x3-x1)*(y2-y1))*((y1-y4)*(x2-x1)+(x4-x1)*(y2-y1))<=0 #intersection inside line P3P4
                    fp2=((y2-y1)*(x4-x3)-(y4-y3)*(x2-x1))*((y4-y3)*(x1-x3)+(y3-y1)*(x4-x3))<=0 #P1, P2 are in the same side of P3P4

                    fp12=fp1*fp2
                    isum[fp12]=isum[fp12]+1
                fp=((isum%2)==0)|((x1==x2)*(y1==y2))
                sind[fp]=0

        #change format
        if fmt==1:
            sindm=np.argmax(sind,axis=1)
            sindm[sind[np.arange(npt),sindm]==0]=-1
            sind=sindm
        elif fmt==0 and npy==1:
            sind=sind[:,0]
    return sind

def find_points_in_polyshp(pt_xy, shapefile_names):
    ind = np.zeros(pt_xy[:, 0].shape)
    for shapefile_name in shapefile_names:
        # shapefile # records = sf.records() # sf.shapeType # len(sf) # s = sf.shape(idx)
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()

        for i, shp in enumerate(shapes):
            poly_xy = np.array(shp.points).T
            print(f'shape {i+1} of {len(shapes)}, {poly_xy[:, 0]}')
            ind += inside_polygon(pt_xy, poly_xy[0], poly_xy[1])  # 1: true; 0: false

    ind = ind.astype('bool')
    return ind


def split_quads(elements=None):  # modified by FY
    '''
    Split quad elements to triangles and append additional elements to element table
    This script can be made much faster by using vector operation instead of the for-loop;
    just append additional elements to the end.
    '''
    from copy import deepcopy

    if elements is None:
        raise Exception('elements should be a numpy array of (np,4)')

    if elements.shape[1] == 3: # already triangles
        return elements
    elif elements.shape[1] != 4:
        raise Exception('elements should be a numpy array of (n,3) or (n,4)')

    triangles = deepcopy(elements)
    quad_idx = ~elements[:, -1].mask
    quads = elements[quad_idx]
    upper_triangle = np.c_[quads[:, 0], quads[:, 1], quads[:, 3], -np.ones((quads.shape[0], 1))]  # last node is masked
    lower_triangle = np.c_[quads[:, 1], quads[:, 2], quads[:, 3], -np.ones((quads.shape[0], 1))]  # fixed a bug where lower triangle is the same as upper triangle

    # replace quads with upper triangle
    triangles[quad_idx, :] = upper_triangle
    # append lower triangle a the end
    triangles = np.ma.concatenate([triangles, lower_triangle], axis=0)
    # mask the last node, because all quads have been changed to triangles
    triangles.mask[:, -1] = True

    # tris = []
    # elements=np.ma.masked_values(elements, -1)  # modified by FY
    # for ele in elements:
    #     ele=ele[~ele.mask]
    #     if len(ele) == 3:
    #         tris.append([ele[0], ele[1], ele[2]])
    #     elif len(ele) == 4:
    #         tris.append([ele[0], ele[1], ele[3]])
    #         tris.append([ele[1], ele[2], ele[3]])

    return triangles[:, :3]  # only return the first 3 nodes of each element

if __name__ == '__main__':
    # Check host and make special arrangement for WCOSS2
    myhost = os.uname()[1]

    # ---------------------------
    my_fillvalue = -99999.0  # used for dry nodes and small disturbance on land/city
    # small_disturbance_ocean_mask_val = -88888.0  # used for small disturbance in ocean
    t0=time()
    # ---------------------------

    # ----------- input arguments ----------------------
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--input_filename', required=True, help='input file in SCHISM format')
    argparser.add_argument('--input_city_identifier_file', help='optional input shapefile defining the urban region, or (on WCOSS2) a txt file containing node ids inside city, where small disturbance will be masked out')
    argparser.add_argument('--datum', default='xGEOID20B', help='Vertical datum of the schism run, defaulting to xGEOID20B.')
    argparser.add_argument('--output_dir', required=True, help='folder holding the outputs of this script.')
    args=argparser.parse_args()
    # ----------- end input arguments ----------------------

    # ----------- process input arguments ----------------------
    input_filename=args.input_filename
    input_city_identifier_file = args.input_city_identifier_file
    datum = args.datum
    output_dir = args.output_dir

    input_fileindex=os.path.basename(input_filename).replace("_", ".").split(".")[1]  # get the file index only

    # decide on the need for a static city mask, sanity checks on the input_city_identifier_file
    if input_city_identifier_file is not None:
        if not os.path.exists(input_city_identifier_file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), input_city_identifier_file)

        if any(tested_host in myhost for tested_host in ["frontera", "viz", "femto", "vortex", "kuro", "gulf"]):
            static_city_mask = False  # search for city mask within polygons of shapefile
            if Path(input_city_identifier_file).suffix != '.shp':
                raise ValueError("When not using static_city_mask, input_city_identifier_file must be a shapefile containing city polygons")
        else:
            static_city_mask = True  # On WCOSS2 or any untested machine, use static mask because it may not have mpl.path
            if Path(input_city_identifier_file).suffix != '.txt':
                raise ValueError("When using static_city_mask, input_city_identifier_file must be a txt file containing node indices inside city")
        print(f'myhost: {myhost}, static_city_mask: {static_city_mask}, city identifier input: {input_city_identifier_file}')
    else:
        print('No city identifier file provided, no city mask will be applied')


    print(f"outputting to {output_dir}")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # ------------------ done processing inputs ------------------

    ds=Dataset(input_filename)
    units=ds['time'].units
    base_date=ds['time'].base_date

    #get coordinates/bathymetry
    x=ds['SCHISM_hgrid_node_x'][:]
    y=ds['SCHISM_hgrid_node_y'][:]
    depth=ds['depth'][:]
    NP=depth.shape[0]

    #get elements and split quads into tris
    elements=ds['SCHISM_hgrid_face_nodes'][:,:]
    tris = split_quads(elements)  # modified by FY
    NE=len(tris)
    print(f'NE is {NE}')
    NV=3

    #get times
    times=ds['time'][:]
    #print(times)
    ntimes=len(times)

    #get elev
    elev=ds['elevation'][:,:]
    elev[elev > 100000] = my_fillvalue

    #mask dry node
    dryFlagNode=ds['dryFlagNode'][:,:]
    elev[dryFlagNode == 1] = my_fillvalue

    #maxelevation
    maxelev=np.max(elev,axis=0)
    idxs=np.argmax(elev,axis=0)
    time_maxelev=times[idxs]

    #find city nodes
    if input_city_identifier_file is None:
        city_node_idx = np.zeros(NP, dtype=bool)  # no city mask
    else:
        if static_city_mask:  # read existing city mask
            city_node_idx = np.loadtxt(input_city_identifier_file, encoding='utf-8').astype(bool)
        else:  # search for city mask within polygons of shapefile, then output to a static file
            output_nodeId_fname = Path(input_city_identifier_file).with_suffix('.node_id.txt')  # save static ids instead of searching (because WCOSS2 has no mpl.path)
            city_node_idx = find_points_in_polyshp(pt_xy=np.c_[x, y], shapefile_names=[input_city_identifier_file])
            np.savetxt(output_nodeId_fname, city_node_idx.astype(int), fmt='%i', encoding='utf-8')

    # disturbance
    maxdist = copy.deepcopy(maxelev)
    land_node_idx = depth < 0
    idxs = np.logical_or(land_node_idx, city_node_idx)
    # cut-off at 0 on normally dry area (depth<0 or city)
    # for maxelev == -99999 (dry), it will also be cut off at 0
    maxdist[idxs] = np.maximum(0, maxelev[idxs] + depth[idxs])

    # set mask for dry nodes
    # this is not necessary because dry nodes are already masked in elev
    # idry=np.zeros(NP)
    # idxs=np.where(maxelev+depth <= 1e-6)
    # print(idxs)
    # maxelev[idxs]=my_fillvalue

    #set mask for small max disturbance (<0.3 m) on land and cities
    small_dist_idx = maxdist < 0.3
    filled_idx=small_dist_idx*(land_node_idx+city_node_idx)
    maxdist[filled_idx]=my_fillvalue
    #set mask for small max disturbance (<0.3 m) in ocean
    # ocean_small_dist_idx=small_dist_idx*(~land_node_idx)
    # maxdist[ocean_small_dist_idx]=small_disturbance_ocean_mask_val

    #get wind speed
    uwind=ds['windSpeedX'][:,:] #,0]
    vwind=ds['windSpeedY'][:,:] #,1]
    idxs=np.where(uwind > 100000)
    uwind[idxs]=my_fillvalue
    vwind[idxs]=my_fillvalue

    ds.close()

    with Dataset(f"{output_dir}/schout_adcirc_{input_fileindex}.nc", "w", format="NETCDF4") as fout:
        #dimensions
        fout.createDimension('time', None)
        fout.createDimension('node', NP)
        fout.createDimension('nele', NE)
        fout.createDimension('nvertex', NV)

        #variables
        fout.createVariable('time', 'f8', ('time',))
        fout['time'].long_name="Time"
#       fout['time'].units = f'seconds since {startdate.year}-{startdate.month}-{startdate.day} 00:00:00 UTC'
#       fout['time'].base_date=f'{startdate.year}-{startdate.month}-{startdate.day} 00:00:00 UTC'
        fout['time'].base_date=base_date
        fout['time'].standard_name="time"
        fout['time'].units=units
        fout['time'][:] = times

        fout.createVariable('x', 'f8', ('node',))
        fout['x'].long_name="node x-coordinate"
        fout['x'].standard_name="longitude"
        fout['x'].units="degrees_east"
        fout['x'].positive="east"
        fout['x'][:]=x

        fout.createVariable('y', 'f8', ('node',))
        fout['y'].long_name="node y-coordinate"
        fout['y'].standard_name="latitude"
        fout['y'].units="degrees_north"
        fout['y'].positive="north"
        fout['y'][:]=y

        fout.createVariable('element', 'i', ('nele','nvertex',))
        fout['element'].long_name="element"
        fout['element'].standard_name="face_node_connectivity"
        fout['element'].start_index=1
        fout['element'].units="nondimensional"
        fout['element'][:]=np.array(tris)

        fout.createVariable('depth', 'f8', ('node',))
        fout['depth'].long_name=f"distance below {datum}"
        fout['depth'].standard_name=f"depth below {datum}"
        fout['depth'].coordinates="y x"  # fixed a bug "time y x"
        fout['depth'].location="node"
        fout['depth'].units="m"
        fout['depth'][:]=depth

        fout.createVariable('zeta_max','f8', ('node',), fill_value=my_fillvalue)
        fout['zeta_max'].standard_name=f"maximum_sea_surface_height_above_{datum}"
        fout['zeta_max'].coordinates="y x"
        fout['zeta_max'].location="node"
        fout['zeta_max'].units="m"
        fout['zeta_max'][:]=maxelev

        fout.createVariable('time_of_zeta_max','f8', ('node',), fill_value=my_fillvalue)
        fout['time_of_zeta_max'].standard_name=f"time_of_maximum_sea_surface_height_above_{datum}"
        fout['time_of_zeta_max'].coordinates="y x"
        fout['time_of_zeta_max'].location="node"
        fout['time_of_zeta_max'].units="sec"
        fout['time_of_zeta_max'][:]=time_maxelev

        fout.createVariable('disturbance_max','f8', ('node',), fill_value=my_fillvalue)
        fout['disturbance_max'].standard_name="maximum_depature_from_initial_condition"
        fout['disturbance_max'].coordinates="y x"
        fout['disturbance_max'].location="node"
        fout['disturbance_max'].units="m"
        fout['disturbance_max'][:]=maxdist

        fout.createVariable('zeta','f8', ('time', 'node',), fill_value=my_fillvalue)
        fout['zeta'].standard_name=f"sea_surface_height_above_{datum}"
        fout['zeta'].coordinates="time y x"
        fout['zeta'].location="node"
        fout['zeta'].units="m"
        fout['zeta'][:,:]=elev

        fout.createVariable('uwind','f8', ('time', 'node',), fill_value=my_fillvalue)
        fout['uwind'].long_name="10m_above_ground/UGRD"
        fout['uwind'].standard_name = "eastward_wind"
        fout['uwind'].coordinates="time y x"
        fout['uwind'].location="node"
        fout['uwind'].units="ms-1"
        fout['uwind'][:,:]=uwind

        fout.createVariable('vwind','f8', ('time', 'node',), fill_value=my_fillvalue)
        fout['vwind'].long_name="10m_above_ground/VGRD"
        fout['vwind'].standard_name = "northward_wind"
        fout['vwind'].coordinates="time y x"
        fout['vwind'].location="node"
        fout['vwind'].units="ms-1"
        fout['vwind'][:,:]=vwind

        fout.title = 'SCHISM Model output'
        fout.source = 'SCHISM model output version v10'
        fout.references = 'http://ccrm.vims.edu/schismweb/'

    print(f'Extraction took {time()-t0} seconds')
