# Data Files
Several of the drivers utilized by the ABC GAN are standard indices based on measurments and can only be found in reference files that are regulaly updated.  These reference files are included in the abcgan-drivers package as package data.  The sources of each of these files is listed below so they can be periodically updated by replacing the files in this directory with updated version.  File names may need to be changed to the default names listed below.

## de421.bsp
Used by `sza` to calculate SZA and ShadHeight and `lunar` to calcualte the moon position and phase.

This file is a JPL ephemeris file that is used internally by skyfield.  This can be updated through the [download function](https://rhodesmill.org/skyfield/api-iokit.html#skyfield.iokit.Loader.download) in the skyfield Loader class.

## dstae.dat
Used by `dst_driver` to generate Dst.

Source URL: https://wdc.kugi.kyoto-u.ac.jp/dstae/index.html

This file can be downloaded from the WDC Kyoto data services website.  Enter the desired date and select "Dst Output" and "IAGA2002-like format", then click "Submit".

## tci_info.txt
Used by `thermo_indx` to generate TCI.

Source URL: https://spaceweather.com/

This file can be downloaded from spaceweather.com.  In the left sidebar under "Thermosphere Climate Index", select the "txt" link from "more data:" to get a complete list of daily TCI to present.

## meiv2.data
Used by `thermo_indx` to generate MEI.

Source URL: https://psl.noaa.gov/enso/mei/

This file can be downloaded from NOAA Physical Sciences Laboratory Multivariate ENSO Index Version 2 (MEI.v2) website.  Click the "Download" link directly above the table.

## MJO_index.dat
Used by `thermo_indx` to generate RMM1 and RMM2.

Source URL: http://www.bom.gov.au/climate/mjo/

This file can be downloaded from the Bureau of Meterology Madden-Julian Oscillation (MJO) website.  Click "RMM Data" under "Download Data".

## lat60N/S90N/S_Means_2hPa-100hPa.h5
Used by `stratospheric_drivers` to calculated drivers related to MERRA2 data.  These files must be generated from original MERRA2 output.

