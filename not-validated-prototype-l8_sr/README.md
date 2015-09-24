## L8SR Version 0.3.1 Release Notes
Release Date: September 23, 2015

### Downloads
L8SR source code

    git clone https://github.com/USGS-EROS/espa-surface-reflectance.git

L8SR auxiliary files

    http://espa.cr.usgs.gov/downloads/auxiliaries/l8sr_auxiliary/l8sr_auxiliary.tar.gz

See git tag [l8_sr-version_0.3.1]

### Installation
  * Install dependent libraries - ESPA product formatter (https://github.com/USGS-EROS/espa-product-formatter)
  * Set up environment variables.  Can create an environment shell file or add the following to your bash shell.  For C shell, use 'setenv VAR "directory"'.
```
    export PREFIX="path_to_directory_for_l8sr_build_data"
```

  * Install baseline auxiliary files and set up the environment variables.
```
    tar -xvzf l8sr_auxiliary.tar.gz
    export L8_AUX_DIR="directory_saved_auxiliary_files"
    (or in c shell use 
    setenv L8_AUX_DIR "directory_saved_auxiliary_files")
```

  * Download (from Github USGS-EROS surface-reflectance project) and install source files. The following build will create a list of executable files under $PREFIX/bin (tested in Linux with the gcc compiler). It will also copy various scripts from the scripts directory up the the $PREFIX/bin directory.
```
    cd l8_sr\c_version\src
    make
    make install
    make clean

    cd l8_sr\landsat_aux\src
    make
    make install
    make clean
```

  * Test - Download Landsat Level 1 files.  Run the do\_l8\_sr Python script in the PREFIX bin directory to run the applications.  Use do\_l8\_sr.py --help for the usage information.  This script requires that your L8SR binaries are in your $PATH or that you have a $BIN environment variable set up to point to the PREFIX bin directory.
```
    convert_lpgs_to_espa --mtl <Landsat_MTL_file> --xml <Landsat_ESPA_XML_file>
    do_l8_sr.py --xml <Landsat_ESPA_XML_file>
```

  * Check output
```
    {scene_name}_toa_*: top-of-atmosphere (TOA) reflectance (or brightness temperature for thermal bands) in internal ESPA file format
    {scene_name}_sr_*: surface reflectance in internal ESPA file format
```

### Dependencies
  * ESPA raw binary and ESPA common libraries from ESPA product formatter and associated dependencies
  * XML2 library
  * Auxiliary data products
    1. LAADS Terra and Aqua CMG and CMA data
    2. CMGDEM HDF file
    3. Various input files and model information provided with the L8SR auxiliary .tar.gz file

### Auxiliary Data Updates
The baseline auxiliary files provided don't include the daily climate data.  In order to generate or update the auxiliary files to the most recent day of year (actually the most current auxiliary files available will be 2-3 days prior to the current day of year do to the latency of the underlying LAADS products) the user will want to run the updatelads.py script available in $PREFIX/bin.  This script can be run with the "--help" argument to print the usage information.  In general the --quarterly argument will reprocess/update all the LAADS data back to 2013.  This is good to do every once in a while to make sure any updates to the LAADS data products are captured.  The --today command-line argument will process the LAADS data for the most recent year.  In general, it is suggested to run the script with --quarterly once a quarter.  Then run the script with --today on a nightly basis.  This should provide an up-to-date version of the auxiliary input data for L8SR.  The easiest way to accomplish this is to set up a nightly and quarterly cron job.

The updatelads script requires a username/password to access the ladssci.nascom.nasa.gov FTP site.  The user will need to contact NASA Contractor/Scientist Sadashiva Devadiga <sadashiva.devadiga-1@nasa.gov> to obtain a username/password for the LAADS FTP site.  In your email explain that you will be using this ftp access to obtain LAADS data for processing Landsat 8 products using the L8SR application provided by the USGS EROS.

updatelads is currently set up to access ESPA_XMLRPC to obtain the ESPA LAADS username/password.  That access will need to be commented out by the user and the user's specific username/password needs to be specified in the script for the username/password.

The following code snippet is how best to handle the updatelads modification.  The __init__ method in updatelads.py should look like the following, and then you will need to put your own LAADS username and password in where it says {put your username/password here}.

```
    def __init__(self):
        # determine the auxiliary directory to store the data
##        xmlrpc = os.environ.get('ESPA_XMLRPC')
##        if xmlrpc is None:
##            msg = "ESPA_XMLRPC environment variable not set... exiting"
##            logger.error(msg)
##            return ERROR
##
##        # get the LAADS username and password
##        try:
##            server = xmlrpclib.ServerProxy(xmlrpc)
##            self.user = server.get_configuration('ladsftp.username')
##            self.password = server.get_configuration('ladsftp.password')
##        except xmlrpclib.ProtocolError, e:
##            msg = "Error connecting to XMLRPC service to fetch credentials: " 
\               
##                "%s" % e
##            logger.error(msg)
##            return ERROR
        self.user = {put your username here}
        self.password = {put your password here}
        print "LAADS FTP username: " + self.user
        print "LAADS FTP password: " + self.password

        # verify that the XMLRPC service returned valid information and
        # the username and password were set in the configuration
##        if len(self.user) <= 0 or len(self.password) <= 0:
##            msg = "Received invalid sized credentials for LAADS FTP from " \
##                "XMLRPC service. Make sure ladsftp.username and " \
##                "ladsftp.password are set in the ESPA_XMLRPC."
##            logger.error(msg)
##            return ERROR
```

### Data Preprocessing
This version of the L8SR application requires the input Landsat products to be in the ESPA internal file format.  After compiling the product formatter raw\_binary libraries and tools, the convert\_lpgs\_to\_espa command-line tool can be used to create the ESPA internal file format for input to the L8SR application.

### Data Postprocessing
After compiling the product-formatter raw\_binary libraries and tools, the convert\_espa\_to\_gtif and convert\_espa\_to\_hdf command-line tools can be used to convert the ESPA internal file format to HDF or GeoTIFF.  Otherwise the data will remain in the ESPA internal file format, which includes each band in the ENVI file format (i.e. raw binary file with associated ENVI header file) and an overall XML metadata file.

### Verification Data

### User Manual

### Product Guide

## Changes From Previous Version
#### Updates on September 23, 2015 - USGS EROS
  1. Fixed a bug accessing the 9x9 window in the land/water mask array. We were accessing invalid memory if the window was on the edges of the scene.
  2. Fixed a bug accessing the CMG arrays for line+1 and sample+1.  We were accessing invalid memory if the scene was at the right or bottom edge of the CMG array.
  3. Modified the update auxiliary files script (updatelads.py) to retry the file download in the event the wget fails.  Cleaned up a few logger issues in this script as well.

