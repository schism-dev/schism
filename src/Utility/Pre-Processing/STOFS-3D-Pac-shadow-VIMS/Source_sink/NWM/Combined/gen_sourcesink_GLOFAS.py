import os
from datetime import datetime, timedelta
import glob
import argparse
import json

import numpy as np
import pandas as pd
from netCDF4 import Dataset
import scipy.interpolate as interp

def is_leap_year(year):
    if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
        return True
    else:
        return False

def read_featureID_file(filename):

    with open(filename) as f:
        lines = f.readlines()
        feature_ids = []
        for line in lines:
            feature_ids.append(line.split('\n')[0])
    return feature_ids

def write_th_file(dataset, timeinterval, fname, issource=True):

    data = []
    for values, interval in zip(dataset, timeinterval):
        if issource:
            data.append(" ".join([f"{interval:G}", *[f'{x: .4f}' for x in values], '\n']))
        else:
            data.append(" ".join([f"{interval:G}", *[f'{-x: .4f}' for x in values], '\n']))

    with open(fname, 'w+') as fid:
        fid.writelines(data)

def write_mth_file(temp, salinity, fname):

    data = []
    #first record
    dt = 0.
    line = [dt]
    [line.append(x) for x in temp]
    [line.append(x) for x in salinity]
    data.append(line)

    #last record
    dt = 1555200.
    line = [dt]
    [line.append(x) for x in temp]
    [line.append(x) for x in salinity]
    data.append(line)
    newset = np.array(data)
    #data.append(" ".join([f"{dt}", *[f'{x: .1f}' for x in salinity], '\n']))

    #with open(fname, 'w+') as fid:
    #    fid.write(str(line))
    np.savetxt(fname, newset, fmt='%.1f')

def get_aggregated_features(nc_feature_id, features):

    aggregated_features = []
    for source_feats in features:
        aggregated_features.extend(list(source_feats))

    in_file=[]
    for feature in aggregated_features:
        idx=np.where(nc_feature_id == int(feature))[0]
        in_file.append(idx.item())

    in_file_2 = []
    sidx = 0
    for source_feats in features:
        eidx = sidx + len(source_feats)
        #in_file_2.append(in_file[sidx:eidx].tolist())
        in_file_2.append(in_file[sidx:eidx])
        sidx = eidx
    return in_file_2

def streamflow_lookup(file, indexes, threshold=-1e-5):
    nc = Dataset(file)
    streamflow = nc["streamflow"][:]
    streamflow[np.where(streamflow < threshold)] = 0.0
    #change masked value to zero
    streamflow[np.where(streamflow.mask)] = 0.0
    data = []
    for indxs in indexes:
        # Note: Dataset already consideres scale factor and offset.
        data.append(np.sum(streamflow[indxs]))
    nc.close()
    return data

def streamflow_lookup_hindcast(file, indexes, threshold=-1e-5, timeind = None):
    ncid = Dataset(file)
    if timeind is None:
       streamflow = ncid["streamflow"][:]  # forecast uses all times
    else:
        streamflow = ncid["streamflow"][timeind,:]  # hindcast doesn't necessarily
    streamflow[np.where(streamflow < threshold)] = 0.0
    #change masked value to zero
    streamflow = streamflow.filled(0.0)
    data = []
    for indxs in indexes:
        # Note: Dataset already consideres scale factor and offset.
        data.append(np.sum(streamflow[:,indxs],axis=1))
    ncid.close()
    return data

def genGLOFASsource(starttime,enddate,glofasfiles,glofasdir,tag =''):
    ''' returns (zero-based) elementids,time, discharge'''
    preproc = np.load(os.path.join(glofasdir,f'glofas_prep{tag}.npz'))
    Erivers=preproc["Erivers"];idxs=preproc["idxs"];ind0N=preproc["ind0N"];ind1N=preproc["ind1N"]; lonN=preproc["lonN"];latN=preproc["latN"];lon=preproc["lon"];lat=preproc["lat"];areaN=preproc["areaN"];final=preproc["final"]
    xctr = preproc["xctr"]
    yctr= preproc["yctr"]
    elnode = preproc["elnode"]
    ne = preproc["ne"]
    uniqErivers=preproc['uniqErivers']
    singles=preproc['singles']
    uniqEriversInd=preproc['uniqEriversInd'] 
    dups=preproc['dups']
    
    # read in spatial subsets of GLOFAS data
    riverCoordsSubset = np.zeros((Erivers.shape[0],2))
    for fname in glofasfiles:
       ncid = Dataset(fname)
       datalat = ncid.variables['latitude'][:]
       datalon = ncid.variables['longitude'][:]
       if np.min(datalon) < 0 and np.max(lonN) >=180: # mismatch between -180-180 data grid and 0-360 grid during prep
           lonN[lonN >180] -=360
       ilat = len(datalat)  -  np.searchsorted(datalat,latN+0.000001,sorter=np.arange(len(datalat)-1, -1, -1)) # lon grid is in descending order, adjust lat value slightly upward to handle potential low precision rounding error 
       ilat[ilat==datalat.shape[0]] -=1 # force ilat index to be within range, out of range lats get excluded by np.isclose check
       ilon = np.searchsorted(datalon,lonN+0.000001)-1
       ilon[ilon<0] == 0 # force ilon index to be within range, otherwise loops around
       g = np.logical_and(np.isclose(datalon[ilon],lonN), np.isclose(datalat[ilat], latN))
       print(f'found {np.sum(g)} locations in area, out of {ind0N.shape}')
       Rlon = ncid.variables['longitude'][ilon[g]]
       Rlat = ncid.variables['latitude'][ilat[g]]
       riverCoordsSubset[g,:] = np.c_[Rlon,Rlat]
       T = ncid.variables['valid_time'][:]
       time_units = ncid.variables['valid_time'].units
       try:
          GLOFASstart = datetime.strptime(time_units,'seconds since %Y-%m-%dT%H:%M:%S')
       except:
          GLOFASstart = datetime.strptime(time_units,'seconds since %Y-%m-%d')
       relmodstart = (startdate - GLOFASstart).total_seconds()
       relmodend = (enddate - GLOFASstart).total_seconds()
       if len(T.shape)>0:
          timeind = np.argwhere(np.logical_and(T>=relmodstart,T<=relmodend)).ravel()
       else: # glofas file with only a single day
          timeind = np.array([0])
          T = np.array([T.data])
       if len(timeind) == 0: # no overlap of data file and requested time
          if T[0]>=relmodend: 
              timeind = np.array([0])
          elif T[-1] < relmodstart:
              timeind = np.array([-1])
       T = T[timeind]-relmodstart
       discharge = np.zeros((timeind.shape[0],ind0N.shape[0])) * np.nan
       nd = ncid.variables['dis24'].shape
       if nd == 4:
         dis24 = ncid.variables['dis24'][timeind,0,:,:]
       else:
         dis24 = ncid.variables['dis24'][timeind,:,:]
       for i in g:
         discharge[:,i] = dis24[:,ilat[i],ilon[i]]
    
    while T[0] > 0:
       T = np.r_[T[0]-86400,T]
       discharge = np.r_[discharge[0,:][np.newaxis,:], discharge]
       print(f'WARNING: mismatch between start of GloFAS data and run starttime, padding start of files with first available GloFAS data.')
    while T[-1] < relmodend-relmodstart:
       T = np.r_[T,T[-1]+86400]
       discharge = np.r_[discharge,discharge[-1,:][np.newaxis,:]]
       print(f'WARNING: mismatch between end of GloFAS data and run endtime, padding end of files with last available GloFAS data.')
    maxdischarge = np.max(discharge,axis=0)
    
    # merge river values from duplicate elements
    dischargeE = np.zeros((discharge.shape[0],uniqErivers.shape[0]))
    dischargeE[:,singles] = discharge[:,uniqEriversInd[singles]]
    for d in dups:
       ind = np.where(Erivers == uniqErivers[d])[0]
       for i in ind:
           dischargeE[:,d] += discharge[:,i]
    
    empty = np.unique(np.argwhere(np.isnan(dischargeE))[:,1])
    if empty.shape[0]>0:
      dischargeE[:,empty] = 0  # fill missing rivers with 0 flow to avoid null values in .th files
      # return sources, T, dischargeE
    return uniqErivers+1, T, dischargeE
 
def read_source_sink(fn):
    d = np.loadtxt(fn,dtype = int)
    source = d[1:d[0]+1] - 1 # one based to zero based, first value is number of locations
    return source

def  genGLOFASfromClima(startdate,enddate,glofasdir,tag=''):
    '''Clima originally generated for full glofas locations, including in US, NWM removed using excludeNWMregfromGloFASClima'''
    d = np.loadtxt(os.path.join(glofasdir,'vsource.th.clim.redis.noNWM')) # need a version where the climatology has not already been redistributed, need to match GloFAS listed points
    s = startdate.timetuple().tm_yday-1
    e = enddate.timetuple().tm_yday
    if e < s: 
        e += 365 +int(is_leap_year(startdate.year))
    T = d[s:e,0] - d[s,0]
    E = read_source_sink(os.path.join(glofasdir,'source_sink.in.clim.redis.noNWM'))
    discharge = d[s:e,1:]
    return E+1,T, discharge

def excludeNWMregfromGloFASClima(glofasdir = '/sciclone/home/cseaton/data10/pacific/GLOFAS/',tag='RUN31l_NWM500km'):
    '''find climatology elements to exclude based on list of regions
    generates static climatology files, does not need to be redone for each run
    expects clima, gis and prep directories inside glofasdir'''
    import shapely as shl
    dis = np.loadtxt(os.path.join(glofasdir,'clima/vsource.th.clim.redis')) # need a version where the climatology has not already been redistributed, need to match GloFAS listed points
    elementids = read_source_sink(os.path.join(glofasdir,'clima/source_sink.in.clim'))
    preproc = np.load(os.path.join(glofasdir,f'prep/glofas_prep{tag}.npz'))
    regfiles=glob.glob(os.path.join(glofasdir,'gis/NWM*reg'))
    xc = preproc["xctr"][elementids]
    yc = preproc["yctr"][elementids]
    ctr = shl.points(np.c_[xc,yc])
    exclude = np.zeros(xc.shape,dtype=bool)
    for f in regfiles:
        reg =  shl.Polygon(np.loadtxt(f,skiprows=3))
        exclude = np.logical_or(reg.contains(ctr), exclude)
    disI = dis[:,np.r_[True,~exclude]]
    if disI[366,0] == 0:
       disI[366:,0] += disI[365,0] + 86400
    write_th_file(list(disI[:,1:]),list(disI[:,0]),os.path.join(glofasdir,'clima/vsource.th.clim.redis.noNWM')) # version excluding NWM regions
    nsink = 0
    nsource = np.sum(~exclude)
    eid_sources2 = elementids[~exclude]+1
    #write source_sink.in
    with open(os.path.join(glofasdir,'clima/source_sink.in.clim.redis.noNWM'), 'w+') as f:
        f.write('{:<d} \n'.format(nsource))
        for eid in eid_sources2:
            f.write('{:<d} \n'.format(int(eid)))
        f.write('\n')
        f.write('{:<d} \n'.format(nsink))


if __name__ == '__main__':
    '''
    Usage: python gen_sourcesink_GLOFAS.py "yyyy-mm-dd" or "yyyy-mm-dd hh:mm:ss" -N ./
    Run this script in oper_3D/NWM/Combined/ directory. Inputs are in the same directory:
        1. sources_{conus, alaska, hawaii}_global.json
        2. glofas_json.npz file specifying elements and glofas locations to extract
        3. glofas data file: glofas-today.nc
        4. NWM data files: nwm.t00z.medium_range.channel_rt_1.f024.{layer}.nc for AK and conus, similar for HI
    If no realtime GloFAS data available, add -C/--clima argument and provide
        3. glofas climatology: vsource.th.clim.redis.noNWM, source_sink.in.clim.redis.noNWM
    can also specify location for 1,2,3,and 4 as arguments.
    '''

    #input paramters 
    argparser = argparse.ArgumentParser()
    argparser.add_argument('date', type=datetime.fromisoformat, help='input file date')
    argparser.add_argument('-g','--glofasfiles',dest='glofasfiles', type=str, default='glofas-today.nc',help='glofas data files, comma separated list to allow composites across separate downloads. Ignored if --clima is used')
    argparser.add_argument('-G','--glofasdir',dest='glofasdir', default='./', type=str, help='location of preproc-glofas.npz')
    argparser.add_argument('-N','--NWMdir',dest='basepath', default='/sciclone/schism10/lcui01/schism20/ICOGS/Pacific/oper_3D/Source_sink/NWM/Combined', type=str, help='location of NWM data files')
    argparser.add_argument('-j','--jsondir',dest='NWMdir', default='./', type=str, help='location of NWM json files')
    argparser.add_argument('-t','--tag',dest='tag', type=str, default='', help='suffix for static npz filename and output image name prefix, defaults to empty string, mainly useful for testing')
    argparser.add_argument('-o','--outtag',dest='outtag', type=str, default='', help='suffix for output files source_sink.in, vsource.th, msource.th, defaults to empty string, mainly useful for testing and for hindcastsi. If not null, startdate appended')
    argparser.add_argument('-d','--days',dest='rnday', type=int, default=3, help='model rndays, defaults to 3 for forecast')
    argparser.add_argument('-H','--hindcast',dest='fcast', action='store_false',default=True, help='hindcast flag, expects NWM file to contain more data that the run period')
    argparser.add_argument('-C','--clima',dest='clima', action='store_true',default=False, help='Glofas clima flag, reads Glofas from year-long .th file and locations from source_sink.in')
    glofasdir = 'glofas/'
    
    args=argparser.parse_args()
    startdate=args.date
    fcast = args.fcast
    clima= args.clima
    glofasdir = args.glofasdir
    basepath = args.basepath
    NWMdir = args.NWMdir
    tag = args.tag
    outtag = args.outtag
    rnday = timedelta(days=args.rnday)
    if len(outtag) > 0:  
       outtag = f"{outtag}-{startdate.strftime('%Y%m%d')}-{rnday.days}days"
       print(outtag)
    glofasfiles = args.glofasfiles.split(',') # glofas download data
    #startdate = datetime(2022, 3, 29, 0)

    #1. region names
    layers = ['conus', 'alaska', 'hawaii']

    localbasepath = './'
    #generate timevector
    timevector1 = np.arange(startdate, rnday-timedelta(days=1), timedelta(days=1)).astype(datetime)
    timevector2 = np.arange(startdate, rnday+timedelta(hours=1), timedelta(hours=1)).astype(datetime)
    enddate = startdate + rnday
    sources_all = {}
    sinks_all = {}
    eid_sources = []
    eid_sinks = []
    times = []
    dates = []

    for layer in layers:
        print(f'layer is {layer}')
        fname_source = os.path.join(NWMdir,f'sources_{layer}_global.json')
        sources_fid = json.load(open(fname_source))

        #add to the final list
        eid_sources.extend(list(sources_fid.keys()))


        #link data
        if fcast:
           files = glob.glob(f'./cached/nwm*.{layer}.nc')
           #remove old files 
           if files is not None:
              print('Remove old files')
              for f in files:
                  try: 
                      os.remove(f)
                  except OSError as e:
                      print("Error: %s : %s" % (f, e.strerror))
           #link new data
           for i, date in enumerate(timevector2):
               if i == 0:
                   date2 = timevector1[0] - timedelta(days=1)
                   if layer == 'conus' or layer == 'alaska':
                       src = f'{basepath}/{date2.strftime("%Y%m%d")}/nwm.t00z.medium_range.channel_rt_1.f024.{layer}.nc'
                   elif layer == 'hawaii':
                       src = f'{basepath}/{date2.strftime("%Y%m%d")}/nwm.t00z.short_range.channel_rt.f02400.{layer}.nc'
               elif i >= 1 and i <= 24:
                   date2 = timevector1[0]
                   if i == 24:
                       it = f'{int(date.hour)+24:03d}'
                   else:
                       it = f'{int(date.hour):03d}'
                   if layer == 'conus' or layer == 'alaska':
                       src = f'{basepath}/{date2.strftime("%Y%m%d")}/nwm.t00z.medium_range.channel_rt_1.f{it}.{layer}.nc'
                   elif layer == 'hawaii':
                       src = f'{basepath}/{date2.strftime("%Y%m%d")}/nwm.t00z.short_range.channel_rt.f{it}00.{layer}.nc'
               else:
                   date2 = timevector1[1]
                   if i >= 25 and i <= 47:
                       it = f'{int(date.hour):03d}'
                   elif i >= 48 and i <= 71:
                       it = f'{int(date.hour)+24:03d}'
                   else:
                       it = f'{int(date.hour)+48:03d}'
                   if layer == 'conus' or layer == 'alaska':
                       src = f'{basepath}/{date2.strftime("%Y%m%d")}/nwm.t00z.medium_range.channel_rt_1.f{it}.{layer}.nc'
                   elif layer == 'hawaii':
                       src = f'{basepath}/{date2.strftime("%Y%m%d")}/nwm.t00z.short_range.channel_rt.f{it}00.{layer}.nc'
               dst = f'{localbasepath}/cached/nwm.t00z.{date.strftime("%Y%m%d%H")}.{layer}.nc'
               os.symlink(src, dst)
           files = glob.glob(f'./cached/nwm*.{layer}.nc')
           files.sort()
           print(f'file 0 is {files[0]}')
           nc_fid0 = Dataset(files[0])["feature_id"][:]
           src_idxs = get_aggregated_features(nc_fid0, sources_fid.values())
   
           sources = []
           for fname in files:
               ds = Dataset(fname)
               ncfeatureid=ds['feature_id'][:]
               if not np.all(ncfeatureid == nc_fid0):
                   print(f'Indexes of feature_id are changed in  {fname}')
                   src_idxs=get_aggregated_features(ncfeatureid, sources_fid.values())
                   nc_fid0 = ncfeatureid
   
               sources.append(streamflow_lookup(fname, src_idxs))
   
               model_time = datetime.strptime(ds.model_output_valid_time, "%Y-%m-%d_%H:%M:%S")
               if layer == 'conus':
                   dates.append(str(model_time))
                   times.append((model_time - startdate).total_seconds())
               ds.close()
           sources_all[layer] = np.array(sources)
        else: # hindcast processing
           files = glob.glob(f'./cached/NWMv3.0_{layer}*.nc')
           #remove old files 
           if files is not None:
              print('Remove old files')
              for f in files:
                  try: 
                      os.remove(f)
                  except OSError as e:
                      print("Error: %s : %s" % (f, e.strerror))
           HIstartdate = datetime(startdate.year - 5,startdate.month, startdate.day)
           HIenddate = datetime(enddate.year - 5,enddate.month, enddate.day)
           for src in glob.glob(f'{basepath}/NWMv3.0_{layer}*-01.nc'):
               dst = f'{localbasepath}/cached/{os.path.basename(src)}'
               os.symlink(src, dst)
           fname = glob.glob(f'./cached/NWMv3.0_{layer}*-01.nc')[0]
           ds = Dataset(fname)
           nc_fid0 = ds["feature_id"][:]
           src_idxs = get_aggregated_features(nc_fid0, sources_fid.values())
   
           ncfeatureid=ds['feature_id'][:]
           if not np.all(ncfeatureid == nc_fid0):
               src_idxs=get_aggregated_features(ncfeatureid, sources_fid.values())
               nc_fid0 = ncfeatureid
           if layer == 'hawaii': # use 2012, data is per minute
              nwmstarttime = datetime.strptime(ds.variables['time'].units, "minutes since %Y-%m-%dT%H:%M:%S")
              timescaler = 60
              relmodend = (HIenddate - nwmstarttime).total_seconds()/timescaler
              relmodstart = (HIstartdate - nwmstarttime).total_seconds() / timescaler
           else:
              nwmstarttime = datetime.strptime(ds.variables['time'].units, "hours since %Y-%m-%dT%H:%M:%S")
              timescaler = 3600
              relmodend = (enddate - nwmstarttime).total_seconds()/timescaler
              relmodstart = (startdate - nwmstarttime).total_seconds() / timescaler
           timeind = np.logical_and(ds.variables['time'][:] >= relmodstart, ds.variables['time'][:] <=relmodend)
           sources = streamflow_lookup_hindcast(fname, src_idxs,timeind = timeind)
           if layer == 'conus':
              times = list((ds.variables['time'][timeind].data - relmodstart)*timescaler)
              dates = [str(startdate + timedelta(seconds = d)) for d in times]
           if layer == 'hawaii': # historical data is 15 minute, subsample
              sources_all[layer] = np.array(sources)[:,::4].T
           else:
              sources_all[layer] = np.array(sources).T
    # Hawaii historical NWM data ends before the end of 2018, pad the remainder of the year with the last value
    if  sources_all['hawaii'].shape[0] < sources_all['conus'].shape[0]:
        a = sources_all['hawaii']
        n = sources_all['conus'].shape[0]
        sources_all['hawaii'] =  np.r_[a,np.repeat(np.mean(a,axis=0)[np.newaxis,:],n-a.shape[0],axis=0)]
    sources = np.concatenate((sources_all['conus'], sources_all['alaska'], sources_all['hawaii']), axis=1) 
 
    #combine with GLOFAS
    if clima:
       Gelementids,Gtime, Gdischarge = genGLOFASfromClima(startdate,enddate,glofasdir,tag=tag)
    else:
       Gelementids,Gtime, Gdischarge = genGLOFASsource(startdate,enddate,glofasfiles,glofasdir,tag=tag)

    match = np.searchsorted(times,Gtime)
    good = match < len(times)
    Gint = interp.interp1d(Gtime,Gdischarge,axis=0)
    Gintdis = Gint(times)
    df_nwm = pd.DataFrame(data=sources, columns=np.array(eid_sources), index=np.array(dates))
    df_glofas = pd.DataFrame(data=Gintdis,columns=np.array(Gelementids), index=np.array(dates))
    df_nwm_transposed = df_nwm.T
    df_glofas_transposed = df_glofas.T

    #concat two dataframe
    df_nwm_transposed.reset_index(inplace=True)
    df_glofas_transposed.reset_index(inplace=True)

    df_source = pd.concat([df_nwm_transposed, df_glofas_transposed])    

    #Combine redundatn elem
    df_source_final = df_source.groupby(by='index', sort=False).sum()
    
    #write to file
    sources2 = df_source_final.T.values
    write_th_file(sources2, times, f'vsource{outtag}.th', issource=True)

    #write msource.th
    eid_sources2 = df_source_final.index.values
    nsource = eid_sources2.shape[0]
    temp = np.full(nsource, -9999.0)
    salt = np.full(nsource, 0.0)
    write_mth_file(temp, salt, f'msource{outtag}.th')

    nsink = 0
    #write source_sink.in
    with open(f'source_sink{outtag}.in', 'w+') as f:
        f.write('{:<d} \n'.format(nsource))
        for eid in eid_sources2:
            f.write('{:<d} \n'.format(int(eid)))
        f.write('\n')

        f.write('{:<d} \n'.format(nsink))
######################################
# model time range from param.nml

