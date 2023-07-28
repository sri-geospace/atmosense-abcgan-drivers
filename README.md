# atmosense-abcgan-drivers
Software to generate the drivers for ABCGAN

## Installation
This package can be installed with pip

```
pip install git+https://github.com/sri-geospace/atmosense-abcgan-drivers.git
```

The following external dependencies are required, some of which may need to be installed manually.
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)
- [h5py](https://www.h5py.org/)
- [flipchem](https://github.com/amisr/flipchem)
- [apexpy](https://apexpy.readthedocs.io/en/latest/)
- [skyfield](https://rhodesmill.org/skyfield/)
- [spacepy](https://spacepy.github.io/)
- [pytz](https://pythonhosted.org/pytz/)

## Usage
This package can be used to generate a full set of drivers both within other python scripts and as a stand-alone command line program.

### Python Script

In the following example, drivers are calculated once an hour on January 1, 2020.

```python
import abcdrivers
import datetime as dt

starttime = dt.datetime(2020,1,1,0,0,0)
endtime = dt.datetime(2020,1,2,0,0,0)
timestep = 1.0  # in decimal hours

site_lat = 65.5
site_lon = -147.5

# generate an array of times to calculate drivers at
times, utime = abcdrivers.generate_time_array(starttime, endtime, timestep)

# calculate drivers for each time
drivers = abcdrivers.collect_drivers(times, site_lat, site_lon)

```

Note that it `collect_drivers()` can take in any array of `datetime.datetime` objects and it is not necessary to use the `generate_time_array()` function if the time array is created elsewhere.

### Command Line

In order to generate driver from the command line, you must create a config file following the format of the example provided in the project root directory.

```
abcgan-drivers config.ini
```

Drivers will be saved to an output hdf5 file along with some other information about the processing.  The output file name is specified in the config file.

## Time Range Limitations

Many drivers are extracted from data files included as package data.  Consequentially, this package can only calculate drivers in the time ranges these data files cover.  To calculate drivers for recent times, the data files may need to be manually updated (see `abcdrivers/data_files/README.md`).  This package cannot calculate drivers for future times.
