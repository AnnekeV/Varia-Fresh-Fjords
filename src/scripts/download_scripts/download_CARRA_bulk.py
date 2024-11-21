### 
# This script bulk downloads 3-hourly CARRA data from the CDS datastore, 
# and saves the data in one GRIB file per variable per month in different subfolders sorted per year.
# The code checks if files already exist in the download folders.
# Make sure the parameters are changed ! And always keep the same naming convention
###

import cdsapi
import itertools
import os

c = cdsapi.Client()

# change this
# base_folder = '/Volumes/IMAU_CARRA/CARRA/3-hourly/' 
base_folder = '/Volumes/imau02/rapid/Anneke/CARRA/3-hourly/'
names = 'CARRA'
domains = 'west_domain'
level_types = 'surface_or_atmosphere'
variables = list(['surface_runoff'])
product_types = 'analysis'
years =  list(range(2010,2024)) 
months = list(range(1,12+1))

# no need to change this
for (year, month, variable) in list(itertools.product(years, months, variables)):
    ofile = ".".join([variable,names,domains,'3h',str(year),str(month),'grib']) # similar naming convention as RACMO
    ofolder = base_folder 
    # make year folder if is does not exist
    if not ofolder:
        os.mkdir(ofolder)
    print('----->' + ofolder + ofile)

    # check if data is not already available
    if os.path.isfile(ofolder + ofile): 
        print('file already exists, moving on...')
        continue

    c.retrieve(
        'reanalysis-carra-single-levels',
        {
            'domain': domains,
            'level_type': level_types,
            'variable': variable,
            'product_type': product_types,
            'time': [
                '00:00', '03:00', '06:00',
                '09:00', '12:00', '15:00',
                '18:00', '21:00',
            ],
            'year': str(year),
            'month': str(month),
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'format': 'grib',
        },
        ofolder + ofile)
