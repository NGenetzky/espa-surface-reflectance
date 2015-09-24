## Ledaps Version 2.3.1 Release Notes
Release Date: September 23, 2015

### Downloads
Ledaps source code

    git clone https://github.com/USGS-EROS/espa-surface-reflectance.git

Ledaps auxiliary files

    http://espa.cr.usgs.gov/downloads/auxiliaries/ledaps_auxiliary/ledaps_aux.1978-2014.tar.gz

See git tag [ledaps-version_2.3.1]

### Installation
  * Install dependent libraries - ESPA product formatter (https://github.com/USGS-EROS/espa-product-formatter)
  * Set up environment variables.  Can create an environment shell file or add the following to your bash shell.  For C shell, use 'setenv VAR "directory"'.
```
    export PREFIX="path_to_directory_for_ledaps_build_data"
```

  * Install baseline auxiliary files.  Please note that the original ozone data has data gaps in 1978 (actual data starts on Nov. 1, 1978), 1979 (partial), 1993 (partial), 1994 (partial), 1995 (complete gap), 1996 (partial), 1997 (partial), 1998 (missing DOY 347+), 2008 (missing Sept. 28/29), ...  You will want to run the updatetoms script (described later) on this baseline set of data.  The NASA LEDAPS group has filled some of these larger data gaps by interpolating the missing data. If the ozone data is missing from the NASA ftp site, then the updatetoms script will not try to update that auxiliary file.
```
    tar -xvzf ledaps_aux.1978-2014.tar.gz
```

  * Setup the environment variables for the auxiliary files
```
    export LEDAPS_AUX_DIR="directory_saved_auxiliary_files"
    (or in c shell use 
    setenv LEDAPS_AUX_DIR "directory_saved_auxiliary_files")
```

  * Download (from Github USGS-EROS surface-reflectance project) and install source files. The following build will create a list of executable files under $PREFIX/bin (tested in Linux with the gcc and gfortran compiler). It will also copy various scripts from the scripts directory up the the $PREFIX/bin directory.
```
    cd ledaps\ledapsSrc\src
    make
    make install
    make clean

    cd ledaps\ledapsAncSrc
    make
    make install
    make clean
```

  * Test - Download Landsat Level 1 files.  Run the do\_ledaps Python script in the LEDAPS bin directory to run the applications.  Use do\_ledaps.py --help for the usage information.  This script requires that your LEDAPS binaries are in your $PATH or that you have a $BIN environment variable set up to point to the LEDAPS bin directory.
```
    convert_lpgs_to_espa --mtl <Landsat_MTL_file> --xml <Landsat_ESPA_XML_file>
    do_ledaps.py --xml <Landsat_ESPA_XML_file>
```

  * Check output
```
    {scene_name}_toa_*: top-of-atmosphere (TOA) reflectance in internal ESPA file format (brightness temperatures are _toa_band6*)
    {scene_name}_sr_*: surface reflectance in internal ESPA file format
```

### Dependencies
  * ESPA raw binary and ESPA common libraries from ESPA product formatter and associated dependencies
  * XML2 library
  * Auxiliary data products
    1. NCEP water vapor data
    2. TOMS ozone data
    3. CMGDEM HDF file

### Auxiliary Data Updates
This baseline auxiliary files provided are good into 2014.  In order to update the auxiliary files to the most recent day of year (actually the most current auxiliary files available will be 2-3 days prior to the current day of year do to the latency of the underlying NCEP and TOMS products) the user will want to run the updatencep.py and updatetoms.py scripts available in $PREFIX/bin.  Both scripts can be run with the "--help" argument to print the usage information for each script.  In general the --quarterly argument will reprocess/update all the NCEP/TOMS data back to 1978.  This is good to do every once in a while to make sure any updates to the NCEP or TOMS data products are captured.  The --today command-line argument will process the NCEP/TOMS data for the most recent year.  In general, it is suggested to run the scripts with --quarterly once a quarter.  Then run the scripts with --today on a nightly basis.  This should provide an up-to-date version of the auxiliary input data for LEDAPS.  The easiest way to accomplish this is to set up a nightly and quarterly cron job.

### Data Preprocessing
This version of the LEDAPS application requires the input Landsat products to be in the ESPA internal file format.  After compiling the espa-common raw\_binary libraries and tools, the convert\_lpgs\_to\_espa command-line tool can be used to create the ESPA internal file format for input to the LEDAPS application.

### Data Postprocessing
After compiling the product-formatter raw\_binary libraries and tools, the convert\_espa\_to\_gtif and convert\_espa\_to\_hdf command-line tools can be used to convert the ESPA internal file format to HDF or GeoTIFF.  Otherwise the data will remain in the ESPA internal file format, which includes each band in the ENVI file format (i.e. raw binary file with associated ENVI header file) and an overall XML metadata file.

### Verification Data

### User Manual

### Product Guide

## Changes From Previous Version
#### Updates on September 23, 2015 - USGS EROS
  * lndsr
    1. lndcal flags saturated pixels however lndsr processes those pixels the same all other valid pixels.  Modifying lndsr to detect saturated pixels, flag them in the output product as such, and keep track of the number of saturated pixels in the lndsr stats output to the screen.
  * ledapsAncSrc
    1. Modified the update auxiliary files scripts (updatencep.py and updatetoms.py) to retry the file download in the event the wget fails.  Also updated convert_ozone.c to handle both OMI and pre-OMI lat/long min/max and step values instead of flagging OMI values as “unexpected”.  Added the scripts to the install section of the Makefile.
