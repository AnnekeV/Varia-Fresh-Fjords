import cdsapi
import numpy as np

c = cdsapi.Client()
pathCARRA = "/Volumes/imau02/rapid/Anneke/CARRA/"

for year in np.arange(2022,2024,2):

    c.retrieve(
        'reanalysis-carra-single-levels',
        {
            'domain': 'west_domain',
            'level_type': 'surface_or_atmosphere',
            'variable': [
                'surface_runoff', #'sea_ice_area_fraction', 
            ],
            'product_type': 'analysis',
            'time': [
                '00:00', 
                '03:00', '06:00',    '09:00', 
                '12:00',
                '15:00',
                '18:00', '21:00',
            ],
            'year': [f'{year}',f'{year+1}'],
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
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
            'format': 'netcdf',
        },
        f'{pathCARRA}_Runoff_{year}_.nc')