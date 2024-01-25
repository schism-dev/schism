from datetime import datetime
import requests
from bs4 import BeautifulSoup
import pandas as pd

date = datetime.utcnow()
print(date.strftime("%Y-%m-%d %H:%M:%S"))
url = 'http://aquatic.pyr.ec.gc.ca/realtimebuoys/default.aspx'

resp = requests.get(url)
soup = BeautifulSoup(resp.text, 'html.parser')
temp = float(soup.find("span", {"id":"MainContent_waterTemp"}).text)

df = pd.DataFrame({'date_utc':[date.strftime("%Y-%m-%d %H:00:00")], 'T': [temp]}) 
df.set_index('date_utc', inplace=True)
df.to_csv('FraserRiver_real-time_waterT.csv', mode='a', header=False)
