import argparse
from datetime import datetime, timedelta
import glob
import os

import numpy as np
import pandas as pd

def get_river_obs(df, datevectors, datevectors2):
    #df = pd.read_csv(fname, sep=',', na_values='')
    df.drop(df.columns[[1, 3, 4, 5]],axis=1, inplace=True)
    #df.rename(columns={df.columns[0]: 'date_local', df.columns[1]: 'flow'}, inplace=True)
    ts = pd.to_datetime(pd.Series(df['date'].values))
    ts2 = ts.dt.tz_localize('US/Eastern', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts2.dt.tz_convert('UTC')
    df.insert(1, 'date_utc', ts3)
    #df.drop('date', inplace=True)
    df.set_index('date_utc', inplace=True)
    df_daily = df.resample('D', closed='left').mean()
    #check missing values
    if df_daily['t'].isnull().sum() > 0:
        df_daily['t'] = df_daily['t'].interpolate()
    
    data = []
    for i, dt in enumerate(datevectors):
        data.append(float(df_daily.loc[dt].values[0]))

    for dt in datevectors2[i+1:]:
        data.append(data[-1])
    data.append(data[-1])
    return data

if __name__ == '__main__':
    '''
    python gen_TEM_1_rescue.py startdate enddate
 
    Assume you have all QC_waterT_*.csv files in the current folder.
    The "enddate" should be the second date in the latest filename. For example:
    the latest filename is "QC_waterT_2023-07-19_2023-07-26.csv", then enddate should be 2023-07-26. The
    script will extend 5 days to the further from enddate.
    python gen_TEM_1_rescue.py 2023-07-01 2023-07-26
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
    files = glob.glob('QC_waterT_*.csv')
    files.sort()

    #check date
    date0 = files[0].split('_')[2]
    date1 = files[-1].split('_')[3].split('.')[0]
    if startDT < datetime.strptime(date0, '%Y-%m-%d'):
        raise ValueError(f'startdate {startDT} is ahead of date {date0} in available files!')
    if endDT > datetime.strptime(date1, '%Y-%m-%d'):
        raise ValueError(f'enddate {endDT} exceeds date {date1} in available files!')
    
    df = pd.concat(map(pd.read_csv, files), ignore_index=True)
    df.drop_duplicates(subset='date', keep='last', inplace=True, ignore_index=True)

    temp = {}
    #get st lawrence river
    temp['SL'] = get_river_obs(df, datevectors, datevectors2)

    rivers = ['SL']
    #write file
    data = []
    for i, date in enumerate(datevectors2):
        line = []
        dt = (date - datevectors[0]).total_seconds()
        print(f'time = {dt}')
        line.append(dt)
        for riv in rivers:
            line.append(temp[riv][i])

        data.append(line)

    newset = np.array(data)
    np.savetxt('TEM_1.th', newset, fmt='%.2f') 
