import argparse
from datetime import datetime, timedelta
import os
import glob

import numpy as np
import pandas as pd

def get_river_obs(df, datevectors, datevectors2):
    #df = pd.read_csv(fname, sep=',', na_values='')
    df.drop(df.columns[[0, 2, 3, 4, 5, 7, 8, 9]],axis=1, inplace=True)
    df.rename(columns={df.columns[0]: 'date_local', df.columns[1]: 'flow'}, inplace=True)
    ts = pd.to_datetime(pd.Series(df['date_local'].values))
    #ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts.dt.tz_convert('UTC')
    df.insert(1, 'date_utc', ts3)
    #df.drop('date_local', inplace=True)
    df.set_index('date_utc', inplace=True)

    df_daily = df.resample('D', closed='left').mean()
    #check missing values
    if df_daily['flow'].isnull().sum() > 0:
        df_daily['flow'] = df_daily['flow'].interpolate()

    data = []
    for i, dt in enumerate(datevectors):
        data.append(round(float(df_daily.loc[dt]['flow']), 3))

    for dt in datevectors2[i+1:]:
        data.append(data[-1])
    return data

if __name__ == '__main__':

    '''
    python gen_fluxth.py startdate enddate
 
    Assume you have all St_Lawrence_river_discharge_*.csv files in the current folder.
    The "enddate" should be the date in the latest filename. For example:
    the latest filename is "St_Lawrence_river_discharge_2023-07-26.csv", 
    then enddate should be 2023-07-26. The
    script will extend 5 days to the further from enddate.
    python gen_fluxth.py 2023-07-01 2023-07-26
    '''
    #input paramters 
    argparser = argparse.ArgumentParser()
    argparser.add_argument('startdate', type=datetime.fromisoformat, help='input startdate')
    argparser.add_argument('enddate', type=datetime.fromisoformat, help='input enddate')
    args=argparser.parse_args()
    startDT=args.startdate
    endDT=args.enddate

    datevectors = pd.date_range(start=startDT.strftime('%Y-%m-%d'), end=endDT.strftime('%Y-%m-%d'))
    #datevector2 - real time in TEM_1.th
    endDT2 = endDT + timedelta(days=5)
    datevectors2 = pd.date_range(start=startDT.strftime('%Y-%m-%d'), end=endDT2.strftime('%Y-%m-%d'))

    #combine csv files
    files = glob.glob('St_Lawrence_river_discharge_*.csv')
    files.sort()

    #check date
    date0 = files[0].split('_')[-1].split('.')[0]
    date1 = datetime.strptime(date0, '%Y-%m-%d') - timedelta(days=3)
    date2 = files[-1].split('_')[-1].split('.')[0]
    if startDT < date1:
        raise ValueError(f'startdate {startDT} is ahead of date {date1} in available files!')
    if endDT > datetime.strptime(date2, '%Y-%m-%d'):
        raise ValueError(f'enddate {endDT} exceeds date {date2} in available files!')
    
    df = pd.concat(map(pd.read_csv, files), ignore_index=True)
    df.drop_duplicates(subset='Date', keep='last', inplace=True, ignore_index=True)

    flow = {}
    #get st Lawrence river
    flow['SL'] = get_river_obs(df, datevectors, datevectors2)


    rivers = ['SL']
    #write file
    data = []
    for i, date in enumerate(datevectors2):
        line = []
        dt = (date - datevectors[0]).total_seconds()
        print(f'time = {dt}')
        line.append(dt)
        for riv in rivers:
            line.append(-flow[riv][i])

        data.append(line)

    newset = np.array(data)
    np.savetxt('flux.th', newset, fmt='%.3f') 
