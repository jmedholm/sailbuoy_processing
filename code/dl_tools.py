import xarray as xr
import pandas as pd
import datetime
from tqdm.notebook import tqdm_notebook as tqdm
import gsw

"""
Module to load and work with SailBuoy data.
"""


def load_data(path):
    
    """
    
    Returns parsed data as 'xarray.Dataset'. 
    
    Walks through each line of the data.txt, dropping lines when the sensor is turning on.
    
    Written by Callum Rollo and Johan Edholm 2022-09-08  
    
    """

    cols = list(pd.read_csv(path, nrows = 1, skiprows=1))
    n = len(cols)
    if any("Message" in s for s in cols):
        n -= 1    # The datalogger ends each line with a ",", which is why we need to skip the last column in each file. For GU's data.txt, the datalogger also adds a Message column, which consists of nothing, messing with the code.
    
    data = pd.read_csv(path,sep=',',header=None,skiprows=1,engine='python',usecols = [i for i in range(n-1)])

    good_samples = 0
    bad_samples = 0
    good_idx = []
    
    for i, msg in enumerate(tqdm(data[0],'Checking lines in input file')):
            if ("Power on" in msg):
                bad_samples += 1
            else:
                good_samples += 1
                good_idx.append(i)

    data = data.loc[good_idx]

    print(f"{good_samples} good measurements found in input file")
    print(f"{bad_samples} bad lines found in input file")
    print('')

    tmp = []
    for i in range(len(list(data.columns))):
        tmp1 = data[i].str.split(' = ',expand=True)
        tmp1 = tmp1.rename(columns={1:tmp1[0][0]})
        tmp.append(tmp1)
    df = pd.concat(tmp,axis=1).drop([0],axis=1)
    df = df.where(df != 'NULL ')


    for i in range(len(df['Time'])):
        df['Time'].iloc[i] = datetime.datetime.strptime(df['Time'].values[i], '%d.%m.%Y %H:%M:%S ').strftime('%Y-%m-%d %H:%M:%S')
    df['Time'] = df['Time'].astype('datetime64[s]')

    for var in df.columns[1:]:
        df[var] = df[var].astype(float)


    df = df.set_index('Time')

    return df.to_xarray()

def fix_standard_attrs(ds):
    
    """
    
    Sets attributes for the standard columns from the datalogger, not any sensors.
    
    Written by Johan Edholm 2022-09-08 (https://github.com/jmedholm)  
    
    """
    
    old_name = ['Time','Lat','Long','TTFF','Count',
             'Commands','TxTries','ONT','DiskUsed',
             'Files','I','V','Temperature']

    standard_name = ['time','latitude','longitude','ttff','count',
                      'commands','tries','ont','disk',
                      'files','current','voltage','temperature']

    units = [' ','degrees_north','degrees_east','s','integer',
             'integer','integer','s','integer',
             'integer','A','V','degrees_c']

    long_name = ['Time in seconds, UTC','Latitude','Longitude','Time to first GPS fix','Transmission counter',
                 'Received commands count','Transmission attempts','Time in acquisition mode','Disk used',
                 'Number of transmitted files','Current consumption','Battery voltage','Payload temperature']

    for i in range(len(old_name)):
        ds[old_name[i]].attrs['standard_name'] = standard_name[i]
        ds[old_name[i]].attrs['long_name'] = long_name[i]
        if old_name[i] != 'Time':
            ds[old_name[i]].attrs['units'] = units[i]
        ds = ds.rename({old_name[i]:standard_name[i]})
        
    return ds

def fix_airmar(ds):
    
    """
    
    Sets attributes for the Airmar output, and corrects error values to NaN.
    
    Written by Johan Edholm 2022-09-08 (https://github.com/jmedholm)  
    
    """
    
    old_name = ['AirmarAirTemp',
                'AirmarWindDirection',
                'AirmarWindSpeed',
                'AirmarWindGust',
                'AirmarHeading',
                'AirmarAirFix']

    var_name = ['air_temp',
                'wind_dir',
                'wind_speed',
                'wind_gust',
                'airmar_heading',
                'airmar_gpsfix']

    standard_name = ['air_temperature',
                     'wind_from_direction',
                     'wind_speed',
                     'wind_speed_of_gust',
                     '',
                     '']

    units = ['degrees_c',
             'degree',
             'm s$^{-1}$',
             'm s$^{-1}$',
             'degree',
             'integer']

    long_name = ['Air temperature',
                 'Wind from',
                 'Wind speed',
                 'Wind gust speed',
                 'Heading from sensor´s internal compass',
                 'GPS fix']
    
    ds[old_name] = ds[old_name].where(ds[old_name] != 999999)

    for i in range(len(old_name)):
        ds[old_name[i]].attrs['standard_name'] = standard_name[i]
        ds[old_name[i]].attrs['long_name'] = long_name[i]
        ds[old_name[i]].attrs['units'] = units[i]
        ds[old_name[i]].attrs['installed_date'] = '2018-10-01'
        ds[old_name[i]].attrs['device_name'] = 'AIRMAR 200WX'
        ds[old_name[i]].attrs['serial_number'] = '60065077'
        ds[old_name[i]].attrs['last_calibrated'] = ''
        ds[old_name[i]].attrs['installed_height'] = '0.7'
        ds[old_name[i]].attrs['firmware'] = ''
        ds[old_name[i]].attrs['vendor_name'] = 'AIRMAR'
        ds[old_name[i]].attrs['model_name'] = '200WX'
        ds[old_name[i]].attrs['model_product_page'] = 'https://www.airmar.com/weather-description.html?id=154' 
        ds[old_name[i]].attrs['note'] = 'If the sensor is not powered on long enough the *airmar_gpsfix* will be 0. This means the sensor has not gained a GPS fix and the data transmitted will be of a lower quality.'
        ds = ds.rename({old_name[i]:var_name[i]})
    
    return ds

def fix_dcps(ds):
    
    """
    
    Sets attributes for the DCPS output, and corrects error values to NaN.
    
    Written by Johan Edholm 2022-09-08 (https://github.com/jmedholm)  
    
    """
    
    old_name = ['DCPSStatus',
                'DCPSOnMin',
                'DCPSSpeed',
                'DCPSDirection']

    var_name = ['dcps_flag',
                'dcps_time',
                'dcps_speed',
                'dcps_dir']

    standard_name = ['dcps_flags',
                     'wind_from_direction',
                     'sea_water_velocity',
                     'sea_water_velocity_to_direction']

    units = ['bitvector',
             'minutes',
             'm s$^{-1}$',
             'degree']

    long_name = ['Status flags',
                 'Acquisition period',
                 'Horizontal current speed',
                 'Current direction']
    
    ds[old_name] = ds[old_name].where(ds[old_name] != 999999)

    for i in range(len(old_name)):
        ds[old_name[i]].attrs['standard_name'] = standard_name[i]
        ds[old_name[i]].attrs['long_name'] = long_name[i]
        ds[old_name[i]].attrs['units'] = units[i]
        ds[old_name[i]].attrs['installed_date'] = '2020-01'
        ds[old_name[i]].attrs['device_name'] = 'DCPS 5400'
        ds[old_name[i]].attrs['serial_number'] = '472'
        ds[old_name[i]].attrs['last_calibrated'] = ''
        ds[old_name[i]].attrs['installed_height'] = '-0.2'
        ds[old_name[i]].attrs['firmware'] = ''
        ds[old_name[i]].attrs['vendor_name'] = 'Aanderaa'
        ds[old_name[i]].attrs['model_name'] = 'DCPS 5400'
        ds[old_name[i]].attrs['model_product_page'] = 'https://www.aanderaa.com/media/pdfs/d411_aanderaa_dcps.pdf' 
        if old_name[i] == 'DCPSStatus':
            ds[old_name[i]].attrs['comment'] = '0 – No errors, Bit 0 - DCPS_Status_Memory_error, Bit 1 - DCPS_Status_File_Open_error, Bit 2 - DCPS_Status_Timeout_error, Bit 3 - DCPS_Status_GPS_error'
        ds = ds.rename({old_name[i]:var_name[i]})
    
    return ds

def fix_aadi(ds):
    
    """
    
    Sets attributes for the AADI output, and corrects error values to NaN.
    
    Written by Johan Edholm 2022-09-12 (https://github.com/jmedholm)  
    
    """
    
    ds['AADI_Salt'] = gsw.SA_from_SP(gsw.SP_from_C(ds['AADI_Cond'],ds['AADI_Temp'],0),ds['longitude'],ds['latitude'],0)
    
    
    old_name = ['AADI_Cond',
                'AADI_Temp',
               'AADI_Salt']

    var_name = ['ssc',
                'sst',
                'sss']

    standard_name = ['sea_water_electrical_conductivity',
                     'sea_water_temperature',
                     'sea_water_absolute_salinity']

    units = ['mS cm$^{-1}$',
             'degrees_c',
             'g kg$^{-1}$']

    long_name = ['Seawater conductivity',
                 'Seawater temperature',
                 'Seawater salinity']
    
    ds['AADI_Cond'] = ds['AADI_Cond'].where((ds['AADI_Cond'] != -1000000))
    ds['AADI_Temp'] = ds['AADI_Temp'].where((ds['AADI_Temp'] != -237.0))
    ds['AADI_Salt'] = ds['AADI_Salt'].where(ds['AADI_Salt'] > 10)
    
    for i in range(len(old_name)):
        ds[old_name[i]].attrs['standard_name'] = standard_name[i]
        ds[old_name[i]].attrs['long_name'] = long_name[i]
        ds[old_name[i]].attrs['units'] = units[i]
        ds[old_name[i]].attrs['installed_date'] = '2018-10-01'
        ds[old_name[i]].attrs['device_name'] = 'AADI Conductivity sensor 4319'
        ds[old_name[i]].attrs['serial_number'] = '1758'
        ds[old_name[i]].attrs['last_calibrated'] = ''
        ds[old_name[i]].attrs['installed_height'] = '-0.7'
        ds[old_name[i]].attrs['firmware'] = ''
        ds[old_name[i]].attrs['vendor_name'] = 'Aanderaa'
        ds[old_name[i]].attrs['model_name'] = 'Conductivity sensor 4319'
        ds[old_name[i]].attrs['model_product_page'] = 'https://www.aanderaa.com/media/pdfs/d369_aanderaa_conductivity_sensor_4319.pdf' 
        ds = ds.rename({old_name[i]:var_name[i]})
    
    return ds