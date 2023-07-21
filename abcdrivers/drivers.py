import datetime as dt
import numpy as np
import h5py
import flipchem
from apexpy import Apex
from skyfield.api import Topos, load
from skyfield.framelib import itrs
from skyfield.almanac import moon_phase
import spacepy.coordinates as spc
import spacepy.time as spt
import pytz
import pkg_resources
from importlib.resources import contents
from importlib.resources import open_binary
from importlib.resources import open_text
from scipy.interpolate import Akima1DInterpolator

def time_convert(time_array, glat, glon):
    # uses apexpy

    # calculate magnetic local time
    A = Apex(time_array[0])
    mlat, mlt = A.convert(glat, glon, 'geo', 'mlt', datetime=np.array(time_array))

    # calculate solar local time
    slt = np.array([(t.hour+(t.minute+t.second/60.)/60.+glon/15.)%24 for t in time_array])

    return {'MLT':mlt, 'SLT':slt, 'MLAT':np.full(mlt.shape,mlat)}


def sza(time_array, glat, glon):
    # uses skyfield

    RE = 6371.2*1000.
    DS = 149.6e6*1000.

    f = pkg_resources.resource_filename(__name__, 'data_files/de421.bsp')
    ephemeris = load(f)
    ts = load.timescale()

    times = ts.from_datetimes([t.replace(tzinfo=pytz.UTC) for t in time_array])

    pfrr = ephemeris['earth'] + Topos(latitude_degrees=glat, longitude_degrees=glon)
    alt, az, dist = pfrr.at(times).observe(ephemeris['sun']).apparent().altaz()
    sza = np.pi/2-alt.radians

    shadow_height = RE*(1./np.cos(sza-np.arcsin(RE/DS*np.sin(sza))-np.arccos(RE/DS))-1.)/1000.
    shadow_height[sza<np.pi/2.] = 0.

    return {'SZA':sza*180./np.pi, 'ShadHeight':shadow_height}

def lunar(time_array):

    f = pkg_resources.resource_filename(__name__, 'data_files/de421.bsp')
    ephemeris = load(f)
    ts = load.timescale()

    times = ts.from_datetimes([t.replace(tzinfo=pytz.UTC) for t in time_array])

    # moon phase
    phase = moon_phase(ephemeris, times).radians
    # moon position
    moon_pos_itrs = ephemeris['earth'].at(times).observe(ephemeris['moon']).apparent().frame_xyz(itrs).km

    # convert to GSE coordinates
    moon_pos_sp = spc.Coords(moon_pos_itrs.T, 'GEO', 'car')
    moon_pos_sp.ticks = spt.Ticktock(time_array, 'UTC')
    moon_pos = moon_pos_sp.convert('GSE', 'car').data.T

    return {'moon_x':moon_pos[0], 'moon_y':moon_pos[1], 'moon_z':moon_pos[2], 'moon_phase':phase*180./np.pi}


def geophys_indx(time_array):

    gpi = {'F10.7':[], 'F10.7avg':[], 'ap':[], 'Ap':[]}
    for t in time_array:
        f107, f107a, ap = flipchem.read_geophys(t)
        gpi['F10.7'].append(f107)
        gpi['F10.7avg'].append(f107a)
        gpi['Ap'].append(ap[0])
        gpi['ap'].append(ap[1])

    for k, v in gpi.items():
        gpi[k] = np.array(v)

    return gpi

def thermo_indx(time_array, verbose=True):
    # TCI
    date_table = []
    tci_table = []

    if verbose:
        print(f"reading file tci_info.txt")
    with pkg_resources.resource_stream(__name__, 'data_files/tci_info.txt') as f:
        for line in f:
            ls = line.split()
            ds = ls[0].split(b'-')
            date_table.append(dt.date(int(float(ds[2])),int(float(ds[0])),int(float(ds[1]))))
            tci_table.append(float(ls[1]))

    date_table = np.array(date_table)
    tci_table = np.array(tci_table)
    nday_table = np.array([(d-date_table[0]).days for d in date_table])
    flags = tci_table<-990.
    tci_table[flags] = np.interp(nday_table[flags], nday_table[~flags], tci_table[~flags])

    # find TCI for all requested times
    idx = [np.flatnonzero(t.date()==date_table)[0] for t in time_array]
    tci = tci_table[idx]


    # MEI
    if verbose:
        print(f"reading file meiv2.data.")
    with pkg_resources.resource_stream(__name__, 'data_files/meiv2.data') as f:
        data = np.loadtxt(f, skiprows=1, max_rows=43)
    utime_table = np.array([[(dt.datetime(int(row[0]),m+1,1)-dt.datetime.utcfromtimestamp(0)).total_seconds() for m,_ in enumerate(row[1:])] for row in data]).flatten()
    mei_table = data[:,1:].flatten()
    flags = mei_table<-990.
    mei_table[flags] = np.interp(utime_table[flags], utime_table[~flags], mei_table[~flags])

    # interpolate MEI to requested times
    utime_array = np.array([(t-dt.datetime.utcfromtimestamp(0)).total_seconds() for t in time_array])
    mei = np.interp(utime_array, utime_table, mei_table)


    # MJO
    if verbose:
        print(f"reading file MJO_index.dat")
    with pkg_resources.resource_stream(__name__, 'data_files/MJO_index.dat') as f:
        data = np.loadtxt(f, skiprows=2, usecols=[0,1,2,3,4])
    date_table = np.array([dt.date(int(d[0]), int(d[1]), int(d[2])) for d in data])
    RMM1_table = data[:,3]
    RMM2_table = data[:,4]
    nday_table = np.array([(d-date_table[0]).days for d in date_table])

    flags = RMM1_table>990.
    RMM1_table[flags] = np.interp(nday_table[flags], nday_table[~flags], RMM2_table[~flags])
    flags = RMM2_table>990.
    RMM2_table[flags] = np.interp(nday_table[flags], nday_table[~flags], RMM2_table[~flags])


    # find MJO indicies for all requested times
    idx = [np.flatnonzero(t.date()==date_table)[0] for t in time_array]
    rmm1 = RMM1_table[idx]
    rmm2 = RMM2_table[idx]

    return {'TCI':tci, 'MEI':mei, 'RMM1':rmm1, 'RMM2':rmm2}

def stratospheric_drivers(time_array, verbose=True):
    hemisphere = "N"
    fname = f'lat60{hemisphere}90{hemisphere}_Means_2hPa-100hPa.h5'
    # accesing the data store in the package folders:
    merra2_N = open_binary('abcdata.drivers.data_files', fname)    
    with h5py.File(merra2_N) as fp:
        levs = fp['/levs'][:]
        UnixTime = fp['/UnixTime'][:]
        Tmean = fp['/Tmean'][:]
        Umean = fp['/Umean'][:]
        
    fT = Akima1DInterpolator(x=UnixTime, y=Tmean, axis = 0)
    fU = Akima1DInterpolator(x=UnixTime, y=Umean, axis = 0)
    req_unixtime = [(x - dt.datetime.utcfromtimestamp(0)).total_seconds() for x in time_array]
    newTmean = fT(req_unixtime)
    newUmean = fU(req_unixtime)
    levels = [2,5,10,30,70,100]
    
    out_dict = {}
    for lev0 in levels:
        lev_i = np.where(levs==lev0)[0]
        out_dict.update({f"T{int(lev0)}":np.squeeze(newTmean[:,lev_i])})
        out_dict.update({f"U{int(lev0)}":np.squeeze(newUmean[:,lev_i])})
    
    return out_dict

def dst_driver(time_array, verbose=True):
    with open_text('abcdata.drivers.data_files', "WWW_dstae03507558.dat") as fp:
        lines = fp.readlines()
        
    curr_line = 0
    dts = list()
    dsts = list()
    while curr_line < len(lines)-1:
        if not lines[curr_line][0].isdigit():
            curr_line += 1
            continue
        sdate, stime, sdoy, dst0 = lines[curr_line].strip().split()
        dts.append(dt.datetime.strptime(sdate+stime,"%Y-%m-%d%H:%M:%S.%f"))
        dsts.append(float(dst0))
        curr_line += 1
    dts = np.array(dts)
    dsts = np.array(dsts)
    
    unixtime = [(x - dt.datetime.utcfromtimestamp(0)).total_seconds() for x in dts]
    dsts_interp = Akima1DInterpolator(x=unixtime, y=dsts)
    req_unixtime = [(x - dt.datetime.utcfromtimestamp(0)).total_seconds() for x in time_array]
    return {'dst':dsts_interp(req_unixtime)}

def collect_drivers(time_array, glat, glon, verbose=True):

    drivers = {}
    if verbose:
        print('converting time.')
    drivers.update(time_convert(time_array, glat, glon))
    if verbose:
        print('Working on sza.')
    drivers.update(sza(time_array, glat, glon))
    if verbose:
        print("Working on the lunar drivers.")
    #drivers.update(lunar(time_array))     # temporarily disable becasue spacepy giving Bus Error 11
    if verbose:
        print('Working on geophysical indices.')
    #drivers.update(geophys_indx(time_array))
    if verbose:
        print('Working on Thermo index.')
    drivers.update(thermo_indx(time_array))
    if verbose:
        print('Working on stratospheric drivers.')
    drivers.update(stratospheric_drivers(time_array))
    if verbose:
        print('Working on dst drivers.')
    drivers.update(dst_driver(time_array))
    
    return drivers



def read_config(config_file):

    # read in config file
    config = configparser.ConfigParser()
    config.read(config_file)
    location_short = config['PARAMS']['SITE_SHORT_NAME']
    location_long  = config['PARAMS']['SITE_DESCRIPTION']
    site_lat = float(config['PARAMS']['SITE_LAT'])
    site_lon = float(config['PARAMS']['SITE_LON'])
    site_info = dict(Site_Short=location_short, Site_Description=location_long, Site_Latitude=site_lat, Site_Longitude=site_lon)
    
    outfile = config['METADATA']['OUTFILE']
    version_number = config['METADATA']['VERSION']
    version_description = config['METADATA']['DESCRIPTION']

    timebinlength = float(config['PARAMS']['TIMEBIN'])  # length of time bins in hours
    # consolidate these with ISO date
    startyear =  int(config['PARAMS']['START_YEAR'])
    startmonth = int(config['PARAMS']['START_MONTH'])
    startday =   int(config['PARAMS']['START_DAY'])
    endyear =    int(config['PARAMS']['END_YEAR'])
    endmonth =   int(config['PARAMS']['END_MONTH'])
    endday =     int(config['PARAMS']['END_DAY'])
    starttime = dt.datetime(startyear,startmonth,startday)
    endtime   = dt.datetime(endyear,endmonth,endday)

    return outfile, starttime, endtime, timebinlength, site_info, version_number, version_description


def generate_time_array(starttime, endtime, timebinlength):
    # change so binlength is not hardcoded
    time_array = np.array([starttime+dt.timedelta(hours=h) for h in np.arange(
        (endtime-starttime).total_seconds()/3600)])
    utime = np.array([(t-dt.datetime.utcfromtimestamp(0)).total_seconds()
                      for t in time_array])
    return time_array


def processing_info(version_number, version_description):
    # Copy/Paste this function from utils here
    proc_info = dict(Version=version_number, Description=version_description)
    return proc_info

# save output file
def save_output(outfile, drivers, utime, site_info=None, proc_info=None):

    print(f"Saving HDF5 file: {outfile}")
    with h5py.File(outfile, 'a') as out:

        # write processing metadata to file
        for key, value in proc_info.items():
            out.create_dataset('ProcessingInfo/{}'.format(key),
                               data=np.array(value.encode('utf-8')))

        # write site info to file
        for key, value in site_info.items():
            out.create_dataset('SiteInfo/{}'.format(key),
                               data=np.array(value.encode('utf-8')))

        # initialize datasets for the first time
        tot_time = out.create_dataset('UnixTime', data=utime, chunks=True,
                                      maxshape=(None,), track_order=True)
        tot_time.attrs.create('full name', 'Unix Time')
        tot_time.attrs.create('description', 'Seconds since January 1, 1970')

        for key in drivers.keys():
            ds = out.create_dataset('Drivers/{}'.format(key),
                                    data=drivers[key], track_order=True,
                                    compression='gzip', compression_opts=1)

    print("done.")

def generate_drivers(config_file):

    outfile, starttime, endtime, timebinlength, site_info, version_number, version_description = read_config(configfile)

    time_array = generate_time_array(starttime, endtime, timebinlength)

    processing_info = processing_info(version_number, version_description)

    # calculate drivers for each time bin
    drivers = collect_drivers(time_array, site_info['Site_Latitude'], site_info['Site_Longitude'])


    save_output(outfile, drivers, time_array, site_info, processing_info)


def main(args=None):
    parser=argparse.ArgumentParser(
        description='''Script for generating a file containing only ABC drivers.''')
    parser.add_argument('config_file', help='Config file.')

    parsed_args = parser.parse_args(args)
    print(f"config_file :{parsed_args.config_file}")
    
    
    generate_drivers(parsed_args.config_file)



#def main():
#    starttime = dt.datetime(1998,1,1)
#    endtime = dt.datetime(2020,12,31)
#
#    time_array = np.array([starttime+dt.timedelta(hours=h) for h in np.arange((endtime-starttime).total_seconds()/3600)])
#    unixtime = np.array([(t-dt.datetime.utcfromtimestamp(0)).total_seconds() for t in time_array])
#    out_drive = collect_drivers(unixtime, 65.5, -147.9)
#
#    print(out_drive)
#
#
#if __name__=='__main__':
#    main()
# with h5py.File('sample_drivers.h5', 'w') as h5:
#     h5.create_dataset('UnixTime', data=utime)
#     h5.create_dataset('MagneticLocalTime', data=mlt)
#     for k, v in gpi_out.items():
#         h5.create_dataset(k, data=v)