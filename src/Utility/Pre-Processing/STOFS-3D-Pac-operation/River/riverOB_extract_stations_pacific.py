import pandas as pd

#USGS stations
stations = ['14211720', '14220500', '14243000', '14105700']
with open("usgs_river.txt",'r') as f:
    lines = f.read().splitlines()

start = []
sites = []
for i, line in enumerate(lines):
    if "Data provided" in line: 
        start.append(i)
        sites.append(line.split(' ')[-1])

for isite in range(len(sites)-1):
    if sites[isite] in stations:
        with open(f'usgs_station_{sites[isite]}.txt', 'w+') as f2:
            for i in range(start[isite], start[isite+1]-1):
                if "WARNING" in lines[i]: break
                f2.write(lines[i])
                f2.write('\n')


#USACE
filename = "BON.Flow-Out.Ave.1Hour.1Hour.txt"

with open("OR00001.xml",'r') as f:
    lines = f.read().splitlines()
info = lines[1]
firsttime = info.split(' ')[2].split('=')[-1].replace('"', '')
lasttime = info.split(' ')[3].split('=')[-1].replace('"', '')

ts = pd.date_range(start=firsttime, end=lasttime, freq='h')
nlen = len(ts)

data = []
for line in lines[2:2+nlen]:
    data.append(float(line.split(' ')[0]))

df = pd.DataFrame({'timestamp': ts, 'flow': data})
df.to_csv(filename, index=False)
