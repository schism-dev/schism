import os
from datetime import datetime, timedelta
import glob
import argparse
import json

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from relocate_source_feeder_lean import SourceSinkIn, relocate_sources

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

if __name__ == '__main__':
    '''
    Usage: python gen_sourcesink.py "yyyy-mm-dd" or "yyyy-mm-dd hh:mm:ss"
    Inputs:
        1. layer name  (e.g., 'conus')
        2. work direcotry "basepath"
        3. old_source_sink_in (from NWM)
        4. folder to save outputs e.g., outdir='./')
        5. relocate map file (relocte_map.txt, generated from relocate_source_feeder.py)
        6. source scaling element id and value
        7. folder where nwm*.nc files saved (e.g., ./combine/)

    '''

    #input paramters 
    argparser = argparse.ArgumentParser()
    argparser.add_argument('date', 
        type=datetime.fromisoformat, 
        help="The date format 'YYYY-MM-DD HH:MM:SS'", 
    )
    args=argparser.parse_args()
    startdate=args.date
    #startdate = datetime(2022, 3, 29, 0)

    #1. region name
    layers = ['conus']

    #2. base path
    basepath = '/sciclone/schism10/lcui01/schism20/ICOGS/ICOGS3D/Forecast/v6/Source_sink/v6'

    #3. old_source_sink_in
    old_source_sink_in = 'source_sink.in.before_relocate'

    #4. folder to save outputs
    outdir = './'

    #5. relocate map file
    relocate_map = np.loadtxt(f'./relocate_map.txt', dtype=int)

    #6. source scaling 
    fname_scale = 'source_scale.txt'

    #7. folder where nwm*.nc saved (
    fdir = 'combine'

    sources_all = {}
    sinks_all = {}
    eid_sources = []
    eid_sinks = []
    times = []
    dates = []

    for layer in layers:
        print(f'layer is {layer}')
        fname_source = f'./sources_{layer}.json'
        fname_sink = f'./sinks_{layer}.json'
        sources_fid = json.load(open(fname_source))
        sinks_fid = json.load(open(fname_sink))

        #add to the final list
        eid_sources.extend(list(sources_fid.keys()))
        eid_sinks.extend(list(sinks_fid.keys()))

        #read nc files
        files = glob.glob(f'./{fdir}/nwm*.{layer}.nc')
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

    sources = sources_all['conus']
    sinks = sinks_all['conus'] 
    
    nsource = np.array(eid_sources).shape[0]
    nsink = np.array(eid_sinks).shape[0]
    print(f'nsource is {nsource}')
    print(f'nsink is {nsink}')

    if os.path.exists('source_sink.in.before_relocate') is False:
        print('Writing source_sink.in.before_relocate file...')
        #write source_sink.in
        with open('source_sink.in.before_relocate', 'w+') as f:
            f.write('{:<d} \n'.format(nsource))
            for eid in eid_sources:
                f.write('{:<d} \n'.format(int(eid)))
            f.write('\n')

            f.write('{:<d} \n'.format(nsink))
            for eid in eid_sinks:
                f.write('{:<d} \n'.format(int(eid)))

    #relocate sources
    df_vsources = relocate_sources(
        old_source_sink_in=SourceSinkIn(filename = old_source_sink_in),
        old_vsource=sources,
        times=np.array(times),
        outdir=outdir,
        relocate_map=relocate_map
    )
    
    #source scaling
    with open(fname_scale) as f:
       total =  f.readline().split(' ')[0]
       print(f'Total sources need to be scaled: {total}!')
       for line in f.read().splitlines():
           scale_idx = line.split(',')[-2].strip()
           scale_value = float(line.split(',')[-1])
           print(scale_idx)
           print(scale_value)
    
           print(f'Pre-scaling is {df_vsources[scale_idx]}')
           df_vsources[scale_idx] = df_vsources[scale_idx]*scale_value
           print(f'Post-scaling is {df_vsources[scale_idx]}')

    #write vsource.th
    df_vsources.to_csv(f'{outdir}/vsource.th', index=False, header=False, sep=' ', float_format='%.2f')
