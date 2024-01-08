from datetime import datetime, timedelta

import numpy as np
import pandas as pd

unit_conv = 0.028316846

def get_usgs_obs(fname=None, datevectors=None):
    df = pd.read_csv(fname, skiprows=26,sep='\t', na_values=' ')
    df.drop([0], inplace=True)
    tz = df['tz_cd'].iloc[1]
    #df.drop(columns=['agency_cd', 'site_no', 'tz_cd', '117345_00060_cd'], inplace=True)
    df.drop(df.columns[[0, 1, 3, 5]], axis=1, inplace=True)
    ts = pd.to_datetime(pd.Series(df['datetime'].values))
    ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts2.dt.tz_convert('UTC')
    df.drop(columns=['datetime'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.insert(1, 'date_utc', ts3)
    df.set_index('date_utc', inplace=True)
    df.rename(columns={df.columns[0]: 'flow'}, inplace=True)

    data = []
    for dt in datevectors:
        data.append(round(float(df.loc[dt]['flow'])*unit_conv, 3))
    data.append(data[-1])
    data.append(data[-1])
    return data

def get_columbia_river_obs(fname, datevectors):
    df = pd.read_csv(fname, sep=',', na_values=' ')
    df= df.rename(columns={df.columns[1]: 'flow'})
    #convert PST to UTC (be careful with DST)
    ts = pd.to_datetime(pd.Series(df['timestamp'].values))
    ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts2.dt.tz_convert('UTC')
    df.insert(1, 'date_utc', ts3)

    df.set_index('date_utc', inplace=True)
    df.drop(columns=['timestamp'], inplace=True)
    data = []
    for dt in datevectors:
        data.append(round(float(df.loc[dt]['flow'])*1000*unit_conv, 3))
    data.append(data[-1])
    data.append(data[-1])
    return data

def get_fraser_river_obs(fname, datevectors):
    df = pd.read_csv(fname, sep=',', na_values='')
    df.drop(df.columns[[0, 2, 3, 4, 5, 7, 8, 9]],axis=1, inplace=True)
    df.rename(columns={df.columns[0]: 'date_local', df.columns[1]: 'flow'}, inplace=True)
    ts = pd.to_datetime(pd.Series(df['date_local'].values))
    #ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts.dt.tz_convert('UTC')
    df.insert(1, 'date_utc', ts3)
    #df.drop('date_local', inplace=True)
    df.set_index('date_utc', inplace=True)
    data = []
    breakpoint()
    for dt in datevectors:
        data.append(round(float(df.loc[dt]['flow']), 3))
    data.append(data[-1])
    data.append(data[-1])
    return data

if __name__ == '__main__':
    date = datetime.now() - timedelta(days=1)
    startDT = datetime(date.year, date.month, date.day)
    endDT = startDT + timedelta(days=1)
    datevectors = pd.date_range(start=startDT.strftime('%Y-%m-%d'), end=endDT.strftime('%Y-%m-%d'))
    endDT2 = startDT + timedelta(days=3)
    datevectors2 = pd.date_range(start=startDT.strftime('%Y-%m-%d'), end=endDT2.strftime('%Y-%m-%d'))

    #orders[Willamette, columbia/3, columbia/3, columbia/3, lewis, Cowlitz, fraser]
    rivers = ['Willamette', 'Columbia', 'Columbia', 'Columbia', 'Lewis', 'Cowlitz', 'Fraser']

    #station ids
    station_ids = {'Willamette': '14211720', 'Lewis': '14220500', 'Cowlitz': '14243000'}
    fname_fraser = 'fraser_river_08MF005_hourly.csv'
    fname_columbia = 'BON.Flow-Out.Ave.1Hour.txt'

    flow = {}
    for key, sid in station_ids.items():
        fname = f'usgs_flow_{sid}.csv'
        flow[key] = get_usgs_obs(fname, datevectors)

    #get fraser river
    flow['Fraser'] = get_fraser_river_obs(fname_fraser, datevectors)

    #get columbia river
    flow['Columbia'] = get_columbia_river_obs(fname_columbia, datevectors)


    #write file
    data = []
    for i, date in enumerate(datevectors2):
        line = []
        dt = (date - datevectors[0]).total_seconds()
        print(f'time = {dt}')
        line.append(dt)
        for riv in rivers:
            if riv == 'Columbia':
                line.append(-round(flow[riv][i]/3, 3))
            else:
                line.append(-flow[riv][i])

        data.append(line)

    newset = np.array(data)
    np.savetxt('flux.th', newset, fmt='%.3f') 
       
    
