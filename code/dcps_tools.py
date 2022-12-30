import pandas as pd
import numpy as np
from collections import Counter
import datetime
import xarray as xr
import gsw
import geomag
import matplotlib.dates as mdates
from tqdm.notebook import tqdm_notebook as tqdm



def parse_gps_str(gps_str):
    
    """
    Helper function to parse GPS sentences
    """
    
    parts = gps_str.split(" ")
    datetime_str = f"{parts[-4]} {parts[-3]}" 
    time = parts[-3]
    lat = parts[-2]
    lon = parts[-1]
    gps_datetime = datetime.datetime.strptime(datetime_str, "%d.%m.%Y %H:%M:%S")
    return gps_datetime, lat, lon

def load_dcps(path,desc):
    
    """
    
    Returns parsed DCPS data as a 'xarray.Dataset'. 
    
    Walks through each line of the DCPS.txt with the following logic:
    
    1. Find a valid open GPS
    2. Read a MEASUREMENT block
    3. Read a valid close GPS
    4. Increment number of measurements (that will define size of arrays) by one
    
    If any of the steps fail, it prints a message and starts again from 1.
    
    Written by Callum Rollo 2022-09-02 (https://callumrollo.github.io/)  
    
    """
    
    adcp_raw = pd.read_csv(path, names=["message"])
    
    good_samples = 0
    no_gps_start, corrupt_sample, no_gps_end = [], [], []

    wait_for = "open"
    all_profile_vars, all_cell_vars = (), ()
    for i, msg in enumerate(tqdm(adcp_raw.message,'Checking lines in input file')):
        if ("opened" in msg) & ('GPS' in msg):
            if wait_for != "open":
                no_gps_end.append(i)
            wait_for="measure"

        if ("MEASUREMENT" in msg) & ('*** Processed data' not in msg):
            if wait_for != "measure":
                no_gps_start.append(i)
                wait_for = "open"
            else:
                measure_msg = msg
                wait_for = "close"
                measure_list = measure_msg.split("\t")[3:]
                keys, vals = measure_list[::2], measure_list[1::2]
                # Find all variables in ADCP measurement cells
                for i, key in enumerate(keys):
                    if "Cell" in key:
                        start_cells = i
                        break
                profile_vars = tuple(set(keys[:i]))
                all_profile_vars = tuple(set(profile_vars + all_profile_vars))
                cell_vars = tuple(set(keys[i:]))
                all_cell_vars = tuple(set(cell_vars + all_cell_vars))
        if ("closed" in msg) & ('GPS' in msg):
            if wait_for != "close":
                corrupt_sample.append(i)
            else:
                good_samples+=1   
            wait_for = "open"

    print(f"{good_samples} good measurements found in input file")
    print(f"{len(corrupt_sample)} corrupt measurements found in input file")
    print(f"{len(no_gps_start) + len(no_gps_end)} measurements without GPS fix found in input file")
    print('')
    num_times = good_samples
    measure_list = measure_msg.split("\t")[3:]
    keys, vals = measure_list[::2], measure_list[1::2]
    
    # Count number of cells
    c = Counter(keys)
    num_cells = c.most_common(1)[0][1]
    cell_ids = np.arange(num_cells)
    print(f"ADCP data comprised of {num_cells} cells")
    print('')
    # Create empty arrays for observations. We use number of measurements and number of cells as dimensions
    profile_dict, cell_dict = {}, {}
    profile_array = np.empty(num_times, dtype=float)
    profile_array[:] = np.nan
    cell_array = np.empty((num_cells, num_times), dtype=float)
    cell_array[:] = np.nan

    # Add all values from measurements
    for var in all_profile_vars:
        profile_dict[var] = profile_array.copy()

    for var in all_cell_vars:
        cell_dict[var] = cell_array.copy()

    # Add GPS times
    time_array = np.empty(num_times, dtype=datetime.datetime)
    time_array[:] = datetime.datetime(1970,1,1)
    profile_dict["time_start"] = time_array.copy()
    profile_dict["lat_start"] = profile_array.copy()
    profile_dict["lon_start"] = profile_array.copy()
    profile_dict["time_end"] = time_array.copy()
    profile_dict["lat_end"] = profile_array.copy()
    profile_dict["lon_end"] = profile_array.copy()

    # Extra items that sometimes appear in data. Bad data flags or something?
    cell_dict['*Horizontal Speed[cm/s]'] = cell_array.copy()
    cell_dict['*Direction[Deg.M]'] = cell_array.copy()
    
    wait_for = "open"
    no_gps_start, corrupt_sample, no_gps_end = [], [], []
    t = 0
    for i, msg in enumerate(tqdm(adcp_raw.message,'Creating dataset')):
        if ("opened" in msg) & ('GPS' in msg):
            if wait_for != "open":
                no_gps_end.append(i)
            else:
                gps_start_msg = msg
                gps_dt, lat, lon = parse_gps_str(gps_start_msg)
                profile_dict["time_start"][t] = gps_dt
                profile_dict["lat_start"][t] = lat
                profile_dict["lon_start"][t] = lon
            wait_for="measure"
        if ("MEASUREMENT" in msg) & ('*** Processed data' not in msg):
            if wait_for != "measure":
                no_gps_start.append(i)
                wait_for = "open"
            else:
                measure_msg = msg
                measure_list = measure_msg.split("\t")[3:]
                keys, vals = measure_list[::2], measure_list[1::2]
                c = -1
                for key, val in zip(keys, vals):
                    if key == "Cell Index":
                        c+=1
                    if key in profile_dict.keys():
                        profile_dict[key][t] = val
                    else:
                        cell_dict[key][c, t] = val
                wait_for = "close"
        if ("closed" in msg) & ('GPS' in msg):
            if wait_for != "close":
                corrupt_sample.append(i)
            else:
                gps_close_msg = msg
                gps_dt, lat, lon = parse_gps_str(gps_close_msg)
                profile_dict["time_end"][t] = gps_dt
                profile_dict["lat_end"][t] = lat
                profile_dict["lon_end"][t] = lon 
                t+=1

            wait_for = "open"
            
    print(f"Loaded {good_samples} good measurements to dataset.")
    print(f"Discarded {len(corrupt_sample)} corrupt measurements.")
    print(f"Discarded {len(no_gps_start) + len(no_gps_end)} measurements without GPS fix.")

    data_vars = {}
    for key, val in profile_dict.items():
        data_vars[key] = (["time"], val)
    for key, val in cell_dict.items():
        data_vars[key] = (["cells", "time"], val)
        
    ds = xr.Dataset(
        data_vars=data_vars,
        coords=dict(
            cells=cell_ids,
            time=profile_dict["time_start"] ,
        ),
        attrs=dict(description=desc),
    )

    return ds

def drop_bad_vars(ds):
    
    """
    
    Takes the dataset, and drop all the variables starting with a '*'. This data is of reduced quality according to Aanderaa.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    
    to_drop = []
    keys = list(ds.keys())
    n = len(keys)

    for i in range(n):
        if keys[i][0] == '*':
            to_drop.append(keys[i])
    print(f"Dropped {len(to_drop)} variables starting with * from the dataset.")
    
    return ds.drop(to_drop)

def add_attrs(ds):
    
    """
    
    Updating the dataset with correct names and attributes (from manual). Almost CF compliant.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    
    
    var_list = ['Heading[Deg.M]',
               'Std Dev Heading[Deg.M]',
               'Pitch[Deg]',
               'Roll[Deg]',
               'Abs Tilt[Deg]',
               'Max Tilt[Deg]',
               'Std Dev Tilt[Deg]',
               'Tilt Direction[Deg]',
               'Charge Voltage Vtx1[V]',
               'Charge Voltage Vtx2[V]',
               'Min Input Voltage[V]',
               'Input Voltage[V]',
               'Input Current[mA]',
               'Memory Used[Bytes]',
               'Record State',
               'Ping Count',
               'Cell Index',
               'Cell State1',
               'Cell State2',
               'Horizontal Speed[cm/s]',
               'Direction[Deg.M]',
               'North Speed[cm/s]',
               'East Speed[cm/s]',
               'Vertical Speed[cm/s]',
               'SP Stdev Horizontal[cm/s]',
               'Strength[dB]',
               'Beam1 Speed[cm/s]',
               'Beam2 Speed[cm/s]',
               'Beam3 Speed[cm/s]',
               'Beam4 Speed[cm/s]',
               'Beam1 Strength[dB]',
               'Beam2 Strength[dB]',
               'Beam3 Strength[dB]',
               'Beam4 Strength[dB]',
               'Beam1 Stdev[cm/s]',
               'Beam2 Stdev[cm/s]',
               'Beam3 Stdev[cm/s]',
               'Beam4 Stdev[cm/s]',
               'Cross Difference[cm/s]',
               'lat_start',
               'lat_end',
               'lon_start',
               'lon_end',
               'time_start',
               'time_end',
               ]
    
    ds = ds[var_list]
    
    standard_names = ['heading',
                      'sd_heading',
                      'pitch',
                      'roll',
                      'abs_tilt',
                      'max_tilt',
                      'sd_tilt',
                      'tilt_direction',
                      'charge_voltage_vtx1',
                      'charge_voltage_vtx2',
                      'min_input_voltage',
                      'input_voltage',
                      'input_current',
                      'memory_used',
                      'record_state',
                      'ping_count',
                      'cell_index',
                      'cell_state1',
                      'cell_state2',
                      'horizontal_speed',
                      'direction',
                      'northward_sea_water_velocity',
                      'eastward_sea_water_velocity',
                      'vertical_speed',
                      'sp_sd_horizontal',
                      'strength',
                      'beam1_speed',
                      'beam2_speed',
                      'beam3_speed',
                      'beam4_speed',
                      'beam1_strength',
                      'beam2_strength',
                      'beam3_strength',
                      'beam4_strength',
                      'beam1_sd',
                      'beam2_sd',
                      'beam3_sd',
                      'beam4_sd',
                      'cross_difference',
                      'lat_start',
                      'lat_end',
                      'lon_start',
                      'lon_end',
                      'time_start',
                      'time_end']
    
    variable_names = ['heading',
                      'sd_heading',
                      'pitch',
                      'roll',
                      'abs_tilt',
                      'max_tilt',
                      'sd_tilt',
                      'tilt_direction',
                      'charge_voltage_vtx1',
                      'charge_voltage_vtx2',
                      'min_input_voltage',
                      'input_voltage',
                      'input_current',
                      'memory_used',
                      'record_state',
                      'ping_count',
                      'cell_index',
                      'cell_state1',
                      'cell_state2',
                      'horizontal_speed',
                      'direction',
                      'vel_north_raw',
                      'vel_east_raw',
                      'vel_up_raw',
                      'sp_sd_horizontal',
                      'strength',
                      'b1_speed',
                      'b2_speed',
                      'b3_speed',
                      'b4_speed',
                      'b1_strength',
                      'b2_strength',
                      'b3_strength',
                      'b4_strength',
                      'b1_sd',
                      'b2_sd',
                      'b3_sd',
                      'b4_sd',
                      'correlation',
                      'lat_start',
                      'lat_end',
                      'lon_start',
                      'lon_end',
                      'time_start',
                      'time_end']

    units = ['degree',
             'degree',
             'degree',
             'degree',
             'degree',
             'degree',
             'degree',
             'degree',
             'V',
             'V',
             'V',
             'V',
             'mA',
             'bytes',
             '',
             '',
             '',
             '',
             '',
             'cm s$^{-1}$',
             'degrees.M',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'dB',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'dB',
             'dB',
             'dB',
             'dB',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'cm s$^{-1}$',
             'degrees',
             'degrees',
             'degrees',
             'degrees',
             '',
             '',
            ]

    long_names = ['Heading',
                  'Standard deviation heading',
                  'Pitch',
                  'Roll',
                  'Absolute tilt',
                  'Maximum tilt',
                  'Standard deviation tilt',
                  'Tilt direction',
                  'Charge voltage Vtx1',
                  'Charge voltage vtx2',
                  'Minimum input voltage',
                  'Input voltage',
                  'Input current',
                  'Memory used',
                  'Record state',
                  'Ping count',
                  'Cell index',
                  'Cell state 1',
                  'Cell state2',
                  'Horizontal speed',
                  'Direction',
                  'North speed uncorrected',
                  'East speed uncorrected',
                  'Vertical speed',
                  'Single ping standard deviation horizontal speed',
                  'Strength',
                  'Beam 1 speed',
                  'Beam 2 speed',
                  'Beam 3 speed',
                  'Beam 4 speed',
                  'Beam 1 strength',
                  'Beam 2 strength',
                  'Beam 3 strength',
                  'Beam 4 strength',
                  'Beam 1 standard deviation',
                  'Beam 2 standard deviation',
                  'Beam 3 standard deviation',
                  'Beam 4 standard deviation',
                  'Cross difference',
                  'Sample start latitude',
                  'Sample end latitude',
                  'Sample start longitude',
                  'Sample end longitude',
                  'Sample start time',
                  'Sample end time',
                 ]

    description = ['Averaged heading from one interval, one heading measurement per ping, vector averaged.',
                   'Standard deviation calculation on all heading values from one interval. Indicates how much the sensor rotates around the vertical axis during a measurement interval.',
                   'Pitch angle, average from one interval, one tilt measurement per ping. Pitch is the rotation angle around the x-axis of the sensor (same axis as Transducer 1 and 3).',
                   'Roll angle, average from one interval, one tilt measurement per ping. Roll is the rotation angle around the y-axis of the sensor (same axis as transducer 4 and 2).',
                   'Angle between sensor plane and horizontal plane. Calculates one value per ping from pitch and roll angles, average of all values.',
                   'Maximum absolute tilt from the interval',
                   'Standard deviation tilt from all values of the absolute tilt in the interval. Indicates if the sensor is moving around with variable tilt during the measurement interval.',
                   'A tilt direction is calculated per ping, this is the average from the interval. Gives the direction where the sensor has its largest tilt, with magnetic north as reference.',
                   'The measured voltage to the capacitor on transmitter electronics to Transducer 1 and 2. It should normally be > 4.8V.',
                   'The measured voltage to the capacitor on transmitter electronics to Transducer 3 and 4. It should normally be > 4.8V.',
                   'The minimum input voltage measured while charging the capacitor bank. It should normally be > 6.0V.',
                   'The voltage measured when not charging while awake, averaged.',
                   'The current measured when not charging while awake, averaged.',
                   'Used heap memory.',
                   'A 32-bit status number which provides quality warnings.',
                   'Number of pings executed, can be lower than configured number of pings to be done.',
                   'An index which gives column and cell number.  A 1xxx is column 1, 2xxx is column 2 and 3xxx is column 3.',
                   'The cell state 1 indicates different conditions like for example if the cell has weak signal or cell is inside blanking zone or illegible zone. If zero, everything is ok.',
                   'The cell state 2 is an addition to cell state 1 and indicates different conditions on beam level like for example beam 1, beam 2, beam 3 or beam 4 is inside blanking zone or illegible zone. If zero, everything is ok.',
                   'Horizontal speed calculated from the 4 beams and compensated for tilt, average of all pings over the interval duration.',
                   'The current direction calculated from the 4 beams combined with compass heading to determine the current direction.',
                   'North speed component, average of all pings over the recording interval, uncorrected.',
                   'East speed component, average of all pings over the recording interval, uncorrected.',
                   'Vertical speed component, average of all pings over the recording interval.',
                   'Single ping standard deviation is the standard deviation calculated from all horizontal speed data in the interval. Divide by the square root of number of pings to get the standard deviation for the set.',
                   'Averaged signal from the four beams. The signal strength is calculated for each beam for each ping and averaged at the end of the interval.',
                   'Averaged speed for beam 1.',
                   'Averaged speed for beam 2.',
                   'Averaged speed for beam 3.',
                   'Averaged speed for beam 4.',
                   'Averaged strength for beam 1.',
                   'Averaged strength for beam 2.',
                   'Averaged strength for beam 3.',
                   'Averaged strength for beam 4.',
                   'Standard deviation calculated from all the beam 1 speeds in the interval.',
                   'Standard deviation calculated from all the beam 2 speeds in the interval.',
                   'Standard deviation calculated from all the beam 3 speeds in the interval.',
                   'Standard deviation calculated from all the beam 4 speeds in the interval.',
                   'Calculated as (Beam 1 Speed + Beam 3 Speed) - (Beam 2 Speed + Beam 4 Speed), where the BeamX speeds are tilt compensated. The Cross Difference is calculated for each ping and averaged at the end of the interval. In a homogeneous water flow this value is close to zero.',
                   'GPS Latitude for the start of measurement.',
                   'GPS Latitude for the end of measurement.',
                   'GPS Longitude for the start of measurement.',
                   'GPS Longitude for the end of measurement.',
                   'GPS Time stamp for the start of measurement.',
                   'GPS Time stamp for the end of measurement.',
                  ]
    
    c = 0
    keys = list(ds.keys()) # Array with all the variables. Some of them have the units in the variable name, and we'd like to move this to the attributes instead.
    for var in keys:
        ds[var].attrs['standard_name'] = standard_names[c]
        ds[var].attrs['long_name'] = long_names[c]
        if 'time' not in var: # It didn't like setting units to a time variable
            ds[var].attrs['units'] = units[c]
        ds[var].attrs['description'] = description[c]
        ds = ds.rename({var:variable_names[c]})
        c += 1        
    print(f"Updated names, units, and description for {len(long_names)} variables")
    
    return ds


def get_bearing2(lat1,lat2,lon1,lon2):
        
    """
    
    Helper function to calculate sample heading.
    
    
    """
    

    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    lon1 = np.deg2rad(lon1)
    lon2 = np.deg2rad(lon2)

    dLon = lon2 - lon1;
    y = np.sin(dLon) * np.cos(lat2);
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) *np.cos(lat2) * np.cos(dLon);
    brng = ((np.rad2deg(np.arctan2(y, x)) + 360) % 360)

    return brng


def calc_mag_dec(ds):
    
    """
    
    Calculates magnetic declination for the given latitudes of the dataset.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    

    mag_dec = np.empty(len(ds['time']), dtype=float)

    for i in range(len(ds.time)):
        mag_dec[i] = geomag.declination(ds.lat_start[i],ds.lon_start[i])    
    
    return mag_dec


def calc_heading_speed(ds):
    
    """
    
    Calculates the heading and speed over ground for the Sailbuoy.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    
    print('Calculating heading and speed of the SB...')
    ds['sample_heading'] = get_bearing2(ds['lat_start'],ds['lat_end'],ds['lon_start'],ds['lon_end'])

    dist = np.empty(len(ds['time']), dtype=float)

    for i in range(len(ds['time'])):
        dist[i] = gsw.distance([ds['lon_start'][i],ds['lon_end'][i]],[ds['lat_start'][i],ds['lat_end'][i]])[0]

    ds['sample_time'] = (ds['time_end']-ds['time_start']).astype('timedelta64[s]')

    ds['sample_vel'] = ds['lat_start']
    ds['sample_vel'].values = dist/ds['sample_time'].values.astype('timedelta64[s]').astype(int)

    ds['sample_u'] = ds['sample_vel'] * np.sin(np.deg2rad(ds['sample_heading']))
    ds['sample_v'] = ds['sample_vel'] * np.cos(np.deg2rad(ds['sample_heading']))
    print('Done!')

    return ds
    
    
def corr_mag_dec(ds):
    
    """
    
    Corrects the currents for magnetic declination.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    
    print('Correcting for magnetic declination...')
    ds['mag_dec'] = (('time'), (calc_mag_dec(ds)))
    ds['c_u_mag'] = ds['horizontal_speed'] * np.sin(np.deg2rad(ds['direction'] + ds['mag_dec'])) / 100
    ds['c_v_mag'] = ds['horizontal_speed'] * np.cos(np.deg2rad(ds['direction'] + ds['mag_dec'])) / 100
    print('Done!')
    return ds

def corr_speed_adcp(ds):
    
    """
    
    Corrects the currents for GPS movement of the Sailbuoy.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    
    print('Correcting for GPS movements...')
    ds['c_u_corr'] = ds['c_u_mag'] + (ds['sample_vel'] * np.sin(np.deg2rad(ds['heading'])))
    ds['c_v_corr'] = ds['c_v_mag'] + (ds['sample_vel'] * np.cos(np.deg2rad(ds['heading'])))
    ds['horizontal_speed'] = (ds['c_u_corr']**2 + ds['c_v_corr']**2)**(1/2)
    print('Done!')

    return ds

def corr_speed(ds):
    
    """
    
    Corrects the currents for GPS movement of the Sailbuoy.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    
    print('Correcting for GPS movements...')
    ds['c_u_corr'] = ds['c_u_mag'] + ds['sample_u']
    ds['c_v_corr'] = ds['c_v_mag'] + ds['sample_v']
    ds['horizontal_speed'] = (ds['c_u_corr']**2 + ds['c_v_corr']**2)**(1/2)
    print('Done!')

    return ds


def add_coords(ds):
     
    """
    
    Adding the start latitude and longitude per sample, as well as depth, as additional coordinates for the dataset.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    
    ds = ds.assign_coords(depth=('cells',np.arange(3,83,2))) # Make sure that we have the correct depth values for the cells. Blanking distance is 3 m.
    ds = ds.assign_coords(latitude=('time', ds['lat_start'].data)).assign_coords(longitude=('time', ds['lon_start'].data))
    print('Added lat, lon, and depth as coordinates.')

    return ds


def update_attrs(ds):
    
    """
    
    Takes the dataset, and updates the newly created variables' attributes.
    
    Written by Johan Edholm, 2022-09-02 
    
    """
    
   
    keys = ['sample_heading',
            'sample_time',
            'sample_vel',
            'sample_u',
            'sample_v',
            'c_u_mag',
            'c_v_mag',
            'c_u_corr',
            'c_v_corr',
            'horizontal_speed']

    standard_names = ['',
                      '',
                      '',
                      '',
                      '',
                      'eastward_sea_water_velocity',
                      'northward_sea_water_velocity',
                      'eastward_sea_water_velocity',
                      'northward_sea_water_velocity',
                      'horizontal_sea_water_velocity']

    unit = ['degrees',
            '',
            'm s$^{-1}$',
            'm s$^{-1}$',
            'm s$^{-1}$',
            'm s$^{-1}$',
            'm s$^{-1}$',
            'm s$^{-1}$',
            'm s$^{-1}$',
            'm s$^{-1}$',]

    long = ['Sample heading',
            'Sample time',
            'Sample velocity',
            'Sample east speed',
            'Sample north speed',
            'East speed magneticly declinated',
            'North speed magneticly declinated',
            'East speed corrected',
            'North speed corrected',
            'Horizontal speed corrected']

    desc = ['Heading of the Sailbuoy during the sample',
            'Time taken to complete the sample',
            'Velocity of the Sailbuoy during the sample',
            'East speed of the Sailbuoy during the sample',
            'North speed of the Sailbuoy during the sample',
            'East speed corrected for magnetic declination',
            'North speed corrected for magnetic declination',
            'East speed corrected for magnetic declination and GPS movement',
            'North speed corrected for magnetic declination and GPS movement',
            'Horizontal speed corrected for magnetic declination and GPS movement',]
    c = 0
    for var in keys:
        if standard_names[c] != '':
            ds[var].attrs['standard_name'] = standard_names[c]
        ds[var].attrs['long_name'] = long[c]
        if 'time' not in var: # It didn't like setting units to a time variable
            ds[var].attrs['units'] = unit[c]
        ds[var].attrs['description'] = desc[c]
        
        c += 1     
    
    print(f"Updated names, units, and description for {len(keys)} variables.")
 
    return ds

