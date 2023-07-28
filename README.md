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

## Drivers

| Key | Full Name | Description |
| --- | --------- | ----------- |
| MLT | Magnetic Local Time | Local time calculated used magnetic coordinate system |
| SLT | Solar Local Time | Local time as determined by the sun's position |
| MLAT | Magnetic Latitude | Apex magnetic latitude |
| SZA | Solar Zenith Angle | Angle between the sun and local zenith |
| ShadHeight | Shadow Height | Height of the Earth's shadow in the atmosphere |
| moon_x | Lunar Position (x) | X component of Lunar position in GSE coordinates |
| moon_y | Lunar Position (y) | Y component of Lunar position in GSE coordinates |
| moon_z | Lunar Position (z) | Z component of Lunar position in GSE coordinates |
| moon_phase | Lunar Phase | Lunar phase in degrees |
| F10.7 | F10.7 Solar Radio Flux | Solar rado flux at 10.7 cm (2800 MHz) in s.f.u (10^-22 W m^-2 Hz^-1) |
| F10.7avg | Average F10.7 Solar Radio Flux | Average solar radio flux at 10.7 cm (2800 MHz) in s.f.u (10^-22 W m^-2 Hz^-1) |
| ap | ap Geomagnetic Index | Three hour equivalent planetary amplitude |
| Ap | Daily Ap Index | Daily equivalent planetary amplitude |
| dst | Disturbance storm-time (Dst) index | Index for geomagnetic storm activity |
| TCI | Thermosphere Climate Index | 60-day running average of global cooling power radiated from the thermosphere |
| MEI | Multivariate ESNO Index | Measure of the El Nino/Southern Oscillation (ESNO)/sea surface temperature |
| RMM1 | First MJO Index | First component of the Madden-Julian Oscillation Indices |
| RMM2 | Second MJO Index | Second component of the Madden-Julian Oscillation Indices |
| T* | Stratospheric Temperature | Stratospheric Temperature at pressure level * |
| U* | Stratospheric Wind | Stratospheric Wind at pressure level * |


## Time Range Limitations

Many drivers are extracted from data files included as package data.  Consequentially, this package can only calculate drivers in the time ranges these data files cover.  To calculate drivers for recent times, the data files may need to be manually updated (see `abcdrivers/data_files/README.md`).  This package cannot calculate drivers for future times.
