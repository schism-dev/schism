import json
import requests
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd

data = []

url = 'https://e-navigation.canada.ca/topics/water-levels/central/temperatures-en?type=recent&location=qc'
resp = requests.get(url)
soup = BeautifulSoup(resp.text, 'html.parser')
scripts=soup.find_all("script")
table = scripts[9].contents[0].split('\r\n')[2].split("= [")[1].split("];")[0]

dicts = table.split('},')
total = len(dicts)
if total <= 1: print('No data available!')
while (total > 1):
    for i, elem in enumerate(dicts):
        if i+1 != total:
            data.append(json.loads(elem+'}'))
        else: 
            data.append(json.loads(elem))
    
    df = pd.DataFrame(data)
    df.set_index('date', inplace=True)
    start_date = df.index[0].split(' ')[0]
    print(start_date)
    end_date = df.index[-1].split(' ')[0]
    print(end_date)
    df.to_csv(f'QC_waterT_{start_date}_{end_date}.csv')
