#!/usr/bin/env python3
from datetime import datetime, timedelta
import argparse

import numpy as np
import pandas as pd

def get_river_discharge(fname, datevectors, datevectors2):
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
    for i, dt in enumerate(datevectors):
        print(f'Getting data for day {i+1}:')
        try:
            data.append(round(float(df.loc[dt]['flow']), 3))
        except:
            if i == 0:
                raise KeyError(f'No discharge data for hindcast {dt}, use old flux.th!')
            else:
                print (f"No discharge data for {dt}, use the previous day's data!")

    for dt in datevectors2[i+1:]:
        data.append(data[-1])
 
    return data

if __name__ == '__main__':

    #input paramters 
    argparser = argparse.ArgumentParser()
    argparser.add_argument('date',
        type=datetime.fromisoformat,
        help="The date format 'YYYY-MM-DD HH:MM:SS'",
    )

    args=argparser.parse_args()
    startdate=args.date
    print(f'startdate is {startdate}')

    #date = datetime.now() - timedelta(days=1)
    #startDT = datetime(date.year, date.month, date.day, 12)
    enddate = startdate + timedelta(days=1)
    datevectors = pd.date_range(start=startdate.strftime('%Y-%m-%d %H:00:00'), end=enddate.strftime('%Y-%m-%d %H:00:00'), tz='UTC')
    enddate2 = startdate + timedelta(days=6)
    datevectors2 = pd.date_range(start=startdate.strftime('%Y-%m-%d %H:00:00'), end=enddate2.strftime('%Y-%m-%d %H:00:00'), tz='UTC')

    rivers = ['St_lawrence']

    #fname = 'QC_02OA016_hourly_hydrometric.csv'
    fname = 'river_st_law_obs.csv'
    flow = {}
    #get st. lawrence river
    flow['St_lawrence'] = get_river_discharge(fname, datevectors, datevectors2)


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
