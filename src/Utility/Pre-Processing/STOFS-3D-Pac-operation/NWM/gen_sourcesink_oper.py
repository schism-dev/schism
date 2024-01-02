import os
from datetime import datetime, timedelta
import glob
import argparse
import json

import numpy as np
import pandas as pd
from netCDF4 import Dataset

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

def write_mth_file(temp, salinity, timeinterval, fname):

    data = []
    for interval in timeinterval:
        data.append(" ".join([f"{interval:G}", *[f'{x: .4f}' for x in temp], '\n']))
    for interval in timeinterval:
        data.append(" ".join([f"{interval:G}", *[f'{x: .4f}' for x in salinity], '\n']))

    with open(fname, 'w+') as fid:
        fid.writelines(data)

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

if __name__ == '__main__':
    '''
    Usage: python extract2asci.py "yyyy-mm-dd" or "yyyy-mm-dd hh:mm:ss"
    Inputs are:
        1. sources_{conus, alaska, hawaii}_global.json
           sinks_{conus, alaska, hawaii}_global.json
        2. climatology file "Dai_Trenberth_climatology_1990-2018_hourly.csv"
        3. work direcotry "basepath"
        4. cached/ (all nwm netcdf files)

    '''

    #input paramters 
    #argparser = argparse.ArgumentParser()
    #argparser.add_argument('date', type=datetime.fromisoformat, help='input file date')
    #args=argparser.parse_args()
    #startdate=args.date
    #startdate = datetime.now().date() - timedelta(days=1)
    date = datetime.now() - timedelta(days=1)
    startdate = datetime(date.year, date.month, date.day)
    print(f'startdate is {startdate}')

    #1. region name
    layers = ['conus', 'alaska', 'hawaii']

    #2.
    clima_filename = 'Dai_Trenberth_climatology_1990-2018_hourly.csv'

    #3. base path
    basepath = '/home1/06923/hyu05/work/oper_3D_Pac/river/NWM'

    #generate timevector
    rnday = timedelta(days=3)
    timevector1 = np.arange(startdate, rnday-timedelta(days=1), timedelta(days=1)).astype(datetime)
    timevector2 = np.arange(startdate, rnday+timedelta(hours=1), timedelta(hours=1)).astype(datetime)

    sources_all = {}
    sinks_all = {}
    eid_sources = []
    eid_sinks = []
    times = []
    dates = []

    for layer in layers:
        print(f'layer is {layer}')
        fname_source = f'./sources_{layer}_global.json'
        fname_sink = f'./sinks_{layer}_global.json'
        sources_fid = json.load(open(fname_source))
        sinks_fid = json.load(open(fname_sink))

        #add to the final list
        eid_sources.extend(list(sources_fid.keys()))
        eid_sinks.extend(list(sinks_fid.keys()))


        #link data
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
            dst = f'{basepath}/cached/nwm.t00z.{date.strftime("%Y%m%d%H")}.{layer}.nc'
            os.symlink(src, dst)

        files = glob.glob(f'./cached/nwm*.{layer}.nc')
        files.sort()
        print(f'file 0 is {files[0]}')
        nc_fid0 = Dataset(files[0])["feature_id"][:]
        src_idxs = get_aggregated_features(nc_fid0, sources_fid.values())
        snk_idxs = get_aggregated_features(nc_fid0, sinks_fid.values())

        sources = []
        sinks = []
        for fname in files:
            ds = Dataset(fname)
            ncfeatureid=ds['feature_id'][:]
            if not np.all(ncfeatureid == nc_fid0):
                print(f'Indexes of feature_id are changed in  {fname}')
                src_idxs=get_aggregated_features(ncfeatureid, sources_fid.values())
                snk_idxs=get_aggregated_features(ncfeatureid, sinks_fid.values())
                nc_fid0 = ncfeatureid

            sources.append(streamflow_lookup(fname, src_idxs))
            sinks.append(streamflow_lookup(fname, snk_idxs))

            model_time = datetime.strptime(ds.model_output_valid_time, "%Y-%m-%d_%H:%M:%S")
            if layer == 'conus':
                dates.append(str(model_time))
                times.append((model_time - startdate).total_seconds())
            ds.close()
        sources_all[layer] = np.array(sources)
        sinks_all[layer] = np.array(sinks)

    sources = np.concatenate((sources_all['conus'], sources_all['alaska'], sources_all['hawaii']), axis=1) 
    sinks = np.concatenate((sinks_all['conus'], sinks_all['alaska'], sinks_all['hawaii']), axis=1) 
    print(sources.shape)
    print(sinks.shape)

    #combine with Dai_Trenberth_climatology
    df = pd.read_csv(clima_filename)
    df.set_index('date', inplace=True)
    mask = (pd.to_datetime(df.index) >= dates[0]) & (pd.to_datetime(df.index) <=dates[-1])
    df_forecast = df[mask]
    df_nwm = pd.DataFrame(data=sources, columns=np.array(eid_sources), index=np.array(dates))
    df_nwm_transposed = df_nwm.T
    df_forecast_transposed = df_forecast.T

    #concat two dataframe
    df_nwm_transposed.reset_index(inplace=True)
    df_forecast_transposed.reset_index(inplace=True)

    df_source = pd.concat([df_nwm_transposed, df_forecast_transposed])    

    #Combine redundatn elem
    df_source_final = df_source.groupby(by='index', sort=False).sum()
    
    #write to file
    sources2 = df_source_final.T.values
    write_th_file(sources2, times, 'vsource.th', issource=True)
    write_th_file(sinks, times, 'vsink.th', issource=False)

    #write msource.th
    eid_sources2 = df_source_final.index.values
    nsource = eid_sources2.shape[0]
    print(f'nsource is {nsource}')
    #temp = np.full(nsource, -9999.0)
    #salt = np.full(nsource, 30.0)
    #write_mth_file(temp, salt, times, 'msource.th')

    nsink = np.array(eid_sinks).shape[0]
    #write source_sink.in
    with open('source_sink.in', 'w+') as f:
        f.write('{:<d} \n'.format(nsource))
        for eid in eid_sources2:
            f.write('{:<d} \n'.format(int(eid)))
        f.write('\n')

        f.write('{:<d} \n'.format(nsink))
        for eid in eid_sinks:
            f.write('{:<d} \n'.format(int(eid)))
