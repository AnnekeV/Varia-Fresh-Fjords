import requests
from bs4 import BeautifulSoup as bs
from urllib.parse import urljoin
import pandas as pd
import csv
from tqdm import tqdm 
import re




res = requests.get("https://en.wikipedia.org/wiki/List_of_fjords_of_Greenland")
soup = bs(res.text, "html.parser")
greenland_fjords = {}
greenland_fjords_coordinates = {}
greenland_fjords_lon, greenland_fjords_lat = {}, {}

for link in soup.find_all("a"):
    url = link.get("href", "")
    absolute_url = urljoin("https://en.wikipedia.org", url)
    if "/wiki/" in url:
        try:
            fjord_name = link.text.strip()
            greenland_fjords[fjord_name] = absolute_url
        except:
            pass

# loop over values not keys
for fjord_name, link in tqdm(greenland_fjords.items()):
    try:
        req = requests.get(link).text
        soup = bs(req, 'lxml')
        latitude = soup.find("span", {"class": "latitude"})
        longitude = soup.find("span", {"class": "longitude"})
        greenland_fjords_coordinates[fjord_name] = (latitude.text, longitude.text)
        greenland_fjords_lat[fjord_name] = latitude.text
        greenland_fjords_lon[fjord_name] = longitude.text
    except:
        pass



df = pd.DataFrame(index=greenland_fjords.keys())
df['url'] = greenland_fjords.values()
dfcoords = pd.DataFrame(index=greenland_fjords_coordinates.keys())
dfcoords['coords'] = greenland_fjords_coordinates.values()
dfcoords['lat'] = greenland_fjords_lat.values()
dfcoords['lon'] = greenland_fjords_lon.values()
# now join the two dataframes based on the index
df = pd.merge(df, dfcoords, left_index=True, right_index=True)



def dms_to_dd(dms):
    '''Converts DMS (Degrees, Minutes, Seconds) coordinates to DD (Decimal Degrees) coordinates
    Args:
    dms: str: DMS coordinates
    Returns:
    dd: float: Decimal Degrees coordinates
    '''

    minutes = 0
    seconds = 0
    try: 
        degrees, minutes, seconds, direction = re.split('[°\'"′″]+', dms)
    except:
        try:
            degrees, minutes, direction = re.split('[°\'"′″]+', dms)
        except:
            degrees, direction = re.split('[°\'"′″]+', dms)
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60)
    if direction in ('S','W'):
        dd*= -1
    return dd

# Apply the function to the 'lat' column
df['lat_dd'] = df['lat'].apply(dms_to_dd)
df['lon_dd'] = df['lon'].apply(dms_to_dd)

df.to_csv('../../../data/processed/greenland_fjords.csv')