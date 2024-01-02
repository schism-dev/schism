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
    df.rename(columns={df.columns[0]: 'temp'}, inplace=True)

    data = []
    for dt in datevectors:
        data.append(float(df.loc[dt]['temp']))
    data.append(data[-1])
    data.append(data[-1])
    return data

def get_columbia_river_obs(fname, datevectors):
    df = pd.read_csv(fname, sep=',', na_values=' ')
    df= df.rename(columns={df.columns[1]: 'T-F'})
    #convert PST to UTC (be careful with DST)
    ts = pd.to_datetime(pd.Series(df['timestamp'].values))
    ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts2.dt.tz_convert('UTC')
    df.insert(1, 'date_utc', ts3)

    df.set_index('date_utc', inplace=True)
    df.drop(columns=['timestamp'], inplace=True)
    data = []
    for dt in datevectors:
        data.append(round((float(df.loc[dt]['T-F']) - 32) * 5 / 9, 1))
    data.append(data[-1])
    data.append(data[-1])
    return data

def get_fraser_river_obs(fname, datevectors):
    df = pd.read_csv(fname, sep=',', na_values='')
    #df.drop(df.columns[[0, 2, 3, 4, 5, 7, 8, 9]],axis=1, inplace=True)
    #df.rename(columns={df.columns[0]: 'date_local', df.columns[1]: 'flow'}, inplace=True)
    df['date_utc']  = pd.to_datetime(pd.Series(df['date_utc'].values))
    #ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    #ts3 = ts.dt.tz_convert('UTC')
    #df.insert(1, 'date_utc', ts3)
    #df.drop('date_local', inplace=True)
    df.set_index('date_utc', inplace=True)
    data = []
    for dt in datevectors:
        data.append(float(df.loc[dt].values[0]))
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
    #No observation at Levis and Cowlitz, use Columbia data
    rivers = ['Willamette', 'Columbia', 'Columbia', 'Columbia', 'Columbia', 'Columbia', 'Fraser']

    #station ids
    station_ids = {'Willamette': '14211720'}
    fname_fraser = 'FraserRiver_real-time_waterT.csv'
    fname_columbia_BON = 'BON.Temp-Water.Inst.1Hour.txt'
    fname_columbia_WRNO = 'WRNO.Temp-Water.Inst.1Hour.txt'

    temp = {}
    for key, sid in station_ids.items():
        fname = f'usgs_temp_{sid}.csv'
        temp[key] = get_usgs_obs(fname, datevectors)

    #get fraser river
    temp['Fraser'] = get_fraser_river_obs(fname_fraser, datevectors)

    #get columbia river
    try:
        print(f'Get temperature data  from Bonneville')
        temp['Columbia'] = get_columbia_river_obs(fname_columbia_BON, datevectors)
    except:
        print(f'No temperature data  at Bonneville, try Warrendale')
        temp['Columbia'] = get_columbia_river_obs(fname_columbia_WRNO, datevectors)

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
       
    
