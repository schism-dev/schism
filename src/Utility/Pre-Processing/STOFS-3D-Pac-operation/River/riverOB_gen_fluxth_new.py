from datetime import datetime, timedelta
import pytz
import argparse

import numpy as np
import pandas as pd

unit_conv = 0.028316846

def get_row_number(filename):
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('agency_cd'):
                skiprows = i
                break
    return skiprows

def get_usgs_obs(df=None, datevectors=None, datevectors2=None):
    #df = pd.read_csv(fname, skiprows=26,sep='\t', na_values=' ')
    df.drop([0], inplace=True)
    tz = df['tz_cd'].iloc[1]
    #df.drop(columns=['agency_cd', 'site_no', 'tz_cd', '117345_00060_cd'], inplace=True)
    df.drop(df.columns[[0, 1, 3]], axis=1, inplace=True)
    ts = pd.to_datetime(pd.Series(df['datetime'].values))
    ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts2.dt.tz_convert('UTC')
    df.drop(columns=['datetime'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.insert(1, 'date_utc', ts3)
    df.set_index('date_utc', inplace=True)
    df.rename(columns={df.columns[0]: 'flow'}, inplace=True)
    df.dropna(axis=0, inplace=True)

    data = []
    for i, dt in enumerate(datevectors):
        print(f'Getting data for day {i+1}:')
        try:
            data.append(round(float(df.loc[dt]['flow'])*unit_conv, 3))
        except:
            if i == 0:
                raise KeyError(f'No discharge data for hindcast {dt}, use old flux.th!')
            else:
                print (f"No discharge data for {dt}, use the previous day's data!")
                data.append(data[-1])

    for dt in datevectors2[i+1:]:
        data.append(data[-1])

    return data

def get_columbia_river_obs(fname, datevectors=None, datevectors2=None):
    df = pd.read_csv(fname, sep=',', na_values=' ')
    df= df.rename(columns={df.columns[1]: 'flow'})
    #convert PST to UTC (be careful with DST)
    ts = pd.to_datetime(pd.Series(df['timestamp'].values))
    #ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts.dt.tz_convert('UTC')
    df.insert(1, 'date_utc', ts3)

    df.set_index('date_utc', inplace=True)
    df.drop(columns=['timestamp'], inplace=True)
    data = []
    for i, dt in enumerate(datevectors):
        print(f'Getting data for day {i+1}:')
        try:
            data.append(round(float(df.loc[dt]['flow'])*unit_conv, 3))
        except:
            if i == 0:
                raise KeyError(f'No discharge data for hindcast {dt}, use old flux.th!')
            else:
               print(f"No discharge data for {dt}, use the previous day's data!")
               data.append(data[-1])

    for dt in datevectors2[i+1:]:
        data.append(data[-1])

    return data

def get_fraser_river_obs(fname, datevectors=None, datevectors2=None):
    df = pd.read_csv(fname, sep=',', na_values='')
    df.drop(df.columns[[0, 2, 3, 4, 5, 7, 8, 9]],axis=1, inplace=True)
    df.rename(columns={df.columns[0]: 'date_local', df.columns[1]: 'flow'}, inplace=True)
    ts = pd.to_datetime(pd.Series(df['date_local'].values))
    #ts2 = ts.dt.tz_localize('US/Pacific', ambiguous='infer', nonexistent='shift_forward')
    ts3 = ts.dt.tz_convert('UTC')
    df.insert(1, 'date_utc', ts3)
    #df.drop('date_local', inplace=True)
    df.set_index('date_utc', inplace=True)
    df.dropna(axis=0, inplace=True)

    data = []
    for i, dt in enumerate(datevectors):
        print(f'Getting data for day {i+1}:')
        try:
            data.append(round(float(df.loc[dt]['flow']), 3))
        except:
            if i == 0:
                raise KeyError(f'No discharge data for hindcast {dt}, use old flux.th!')
            else:
                print (f"No discharge data for {dt}, use the previous day's data!")
                data.append(data[-1])
            
    for dt in datevectors2[i+1:]:
        data.append(data[-1])

    return data

if __name__ == '__main__':

    '''
    usage: python gen_fluxth.py YYYY-mm-dd-hh (e.g., 2023-06-20-12)

    '''

    argparser = argparse.ArgumentParser()
    argparser.add_argument('date',
        type=datetime.fromisoformat,
        help="The date format 'YYYY-MM-DD HH:MM:SS'",
    )

    args=argparser.parse_args()
    startDT=args.date
    print(f'startdate is {startDT}')

    endDT = startDT + timedelta(days=1)
    datevectors = pd.date_range(start=startDT.strftime('%Y-%m-%d %H:00:00'), end=endDT.strftime('%Y-%m-%d %H:00:00'), tz=pytz.UTC)
    endDT2 = startDT + timedelta(days=4)
    datevectors2 = pd.date_range(start=startDT.strftime('%Y-%m-%d %H:00:00'), end=endDT2.strftime('%Y-%m-%d %H:00:00'), tz=pytz.UTC)

    #station_ids = {'Willamette': '14211720', 'Lewis': '14220500', 'Cowlitz': '14243000', 'Dalles': '14105700'}

    #orders[Willamette, columbia/3, columbia/3, columbia/3, lewis, Cowlitz, fraser]
    rivers = ['Willamette', 'Columbia', 'Columbia', 'Columbia', 'Lewis', 'Cowlitz', 'Fraser']

    flow = {}
    #Willamette
    print(f'Getting data for Wilamette river...')
    filename = 'usgs_station_14211720.txt'
    rows = get_row_number(filename)
    print(f'skiprows is {rows}')
    df1 = pd.read_csv(filename, skiprows=rows, sep='\t', na_values=' ')
    df1.drop(df1.iloc[:,5:], axis=1, inplace=True)
    flow['Willamette'] = get_usgs_obs(df1, datevectors, datevectors2)

    #read BovillE station for columbia river
    print(f'Getting data for Columbia river...')
    fname_columbia = 'BON.Flow-Out.Ave.1Hour.1Hour.txt'
    flow['Columbia'] = get_columbia_river_obs(fname_columbia, datevectors, datevectors2)

    #Lewis
    print(f'Getting data for Lewis river...')
    filename = 'usgs_station_14220500.txt'
    rows = get_row_number(filename)
    print(f'skiprows is {rows}')
    df3 = pd.read_csv(filename, skiprows=rows, sep='\t', na_values=' ')
    df3.drop(df3.iloc[:,5:], axis=1, inplace=True)
    flow['Lewis'] = get_usgs_obs(df3, datevectors, datevectors2)

    #Cowlitz
    print(f'Getting data for Cowlitz river...')
    filename = 'usgs_station_14243000.txt'
    rows = get_row_number(filename)
    print(f'skiprows is {rows}')
    df4 = pd.read_csv(filename, skiprows=rows, sep='\t', na_values=' ')
    df4.drop(df4.columns[[4, 5, 7, 8, 9, 10, 11, 12, 13]], axis=1, inplace=True)
    flow['Cowlitz'] = get_usgs_obs(df4, datevectors, datevectors2)

    #Fraser river
    print(f'Getting data for Fraser river...')
    fname_fraser = 'BC_08MF005_hourly_hydrometric.csv'
    #get fraser river
    flow['Fraser'] = get_fraser_river_obs(fname_fraser, datevectors, datevectors2)


    #write file
    data = []
    for i, date in enumerate(datevectors2[:-1]):
        line = []
        dt = (date - datevectors[0]).total_seconds()
        print(f'time = {dt}')
        line.append(dt)
        for riv in rivers:
            print(riv)
            if riv == 'Columbia':
                line.append(-round(flow[riv][i]/3, 3))
            else:
                line.append(-flow[riv][i])

        data.append(line)

    newset = np.array(data)
    np.savetxt('flux.th', newset, fmt='%.3f')

