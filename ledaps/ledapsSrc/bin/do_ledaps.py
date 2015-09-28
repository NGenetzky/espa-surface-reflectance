#! /usr/bin/env python
import sys
import os
import re
import commands
import datetime
from optparse import OptionParser

ERROR = 1
SUCCESS = 0


############################################################################
# Description: isLeapYear will determine if the specified year is a leap
# year.
#
# Inputs:
#   year - year to determine if it is a leap year (integer)
#
# Returns:
#     1 - yes, this is a leap year
#     0 - no, this is not a leap year
#
# Notes:
############################################################################
def isLeapYear(year):
    if (year % 4) == 0:
        if (year % 100) == 0:
            if (year % 400) == 0:
                return True
            else:
                return False
        else:
            return True
    else:
        return False


############################################################################
# Description: logIt logs the information to the logfile (if valid) or to
# stdout if the logfile is None.
#
# Inputs:
#   msg - message to be printed/logged
#   log_handler - log file handler; if None then print to stdout
#
# Returns: nothing
#
# Notes:
############################################################################
def logIt(msg, log_handler):
    if log_handler is None:
        print msg
    else:
        log_handler.write(msg + '\n')


#############################################################################
# Created on November 13, 2012 by Gail Schmidt, USGS/EROS
# Created Python script so that lnd* application status values can be checked
# for successful completion and failures can be flagged.
#
# History:
#   Updated on 11/29/2012 by Gail Schmidt, USGS/EROS
#   Created the Ledaps class and added modules to run the LEDAPS code and
#       determine if the proper ancillary data exists for a requested year
#       and year/doy.
#   Updated on 12/6/2012 by Gail Schmidt, USGS/EROS
#   Modified the run_ledaps method to change to the XML directory
#       before running the LEDAPS software, and then return to the current
#       directory upon successful or failure exit.
#   Updated on 3/7/2014 by Gail Schmidt, USGS/EROS
#   Modified the run_ledaps script to work with the input XML file as
#       part of the switch to utilize the ESPA internal raw binary file
#   Updated on 10/09/2014 by Ron Dilley, USGS/EROS
#   Modified to be pep8 compliant.
#
# Usage: do_ledaps.py --help prints the help message
############################################################################
class Ledaps():

    def __init__(self):
        pass

    ########################################################################
    # Description: findAncillary will parse the ANC_PATH directory and verify
    # that the required NCEP REANALYSIS and EP/TOMS ancillary products are
    # available for the specified year and DOY.  If the DOY is not specified,
    # then the entire year will be processed.
    #
    # Inputs:
    #   year - which year to look for ancillary data (integer)
    #   doy - which DOY to look for ancillary data (integer)
    #         (default - parse all days in the specified year)
    #
    # Returns:
    #     None - error occurred while processing the specified year/doy
    #     list of True/False - current value is True if the ancillary data is
    #         available for the current year/doy, False if the ancillary data
    #         does not exist;  if the DOY was specified, then the list will
    #         only have one value
    #
    # Notes:
    #     ANC_PATH points to the base LEDAPS ancillary directory which
    #         contains the REANALYSIS and EP/TOMS subdirectories.
    #######################################################################
    def findAncillary(self, year, doy=-99):
        # determine the ancillary directory to store the data
        ancdir = os.environ.get('ANC_PATH')
        if ancdir is None:
            print "ANC_PATH environment variable not set... exiting"
            return None

        # initialize the doyList to empty and the number of days to 1
        doyList = []
        ndays = 1

        # if doy wasn't specified, then determine the number of days in the
        # specified year. if the specified year is the current year, only
        # process up through today otherwise process through all the days
        # in the year
        if doy == -99:
            now = datetime.datetime.now()
            if year == now.year:
                ndays = now.timetuple().tm_yday
            else:
                if isLeapYear(year) is True:
                    ndays = 366
                else:
                    ndays = 365

        # loop through the number of days and determine if the required
        # ancillary data exists
        for currdoy in range(1, ndays+1):
            # if the DOY was specified then use that value otherwise use
            # the current DOY based on the number of days to be processed
            if doy != -99:
                currdoy = doy

            # pad the DOY with 0s if needed
            if currdoy < 10:
                dayofyear = '00' + str(currdoy)
            elif 9 < currdoy < 100:
                dayofyear = '0' + str(currdoy)
            else:
                dayofyear = str(currdoy)

            # NCEP REANALYSIS file
            ncepFile = ("%s/REANALYSIS/RE_%d/REANALYSIS_%d%s.hdf"
                        % (ancdir, year, year, dayofyear))

            # EP/TOMS file
            tomsFile = ("%s/EP_TOMS/ozone_%d/TOMS_%d%s.hdf"
                        % (ancdir, year, year, dayofyear))
            if os.path.isfile(ncepFile) and os.path.isfile(tomsFile):
                doyList.append(True)
            else:
                doyList.append(False)

        # return the True/False list
        return doyList

    ########################################################################
    # Description: runLedaps will use the parameters passed for xmlfile,
    # logfile, and usebin.  If xmlfile is None (i.e. not specified) then
    # the command-line parameters will be parsed for this information.
    # The LEDAPS applications are then executed on the specified XML file.
    # If a log file was specified, then the output from each LEDAPS application
    # will be logged to that file.
    #
    # Inputs:
    #   xmlfile - name of the Landsat XML file to be processed
    #   process_sr - specifies whether the surface reflectance processing,
    #       should be completed.  True or False.  Default is True, otherwise
    #       the processing will halt after the TOA reflectance products are
    #       complete.
    #   logfile - name of the logfile for logging information; if None then
    #       the output will be written to stdout
    #   usebin - this specifies if the LEDAPS exes reside in the $BIN
    #       directory; if None then the LEDAPS exes are expected to be in
    #       the PATH
    #
    # Returns:
    #     ERROR - error running the LEDAPS applications
    #     SUCCESS - successful processing
    #
    # Notes:
    #   1. The script obtains the path of the XML file and changes
    #      directory to that path for running the LEDAPS code.  If the
    #      xmlfile directory is not writable, then this script exits with
    #      an error.
    #######################################################################
    def runLedaps(self, xmlfile=None, process_sr="True", logfile=None,
                  usebin=None):
        # if no parameters were passed then get the info from the
        # command line
        if xmlfile is None:
            # get the command line argument for the XML file
            parser = OptionParser()
            parser.add_option("-f", "--xml",
                              type="string", dest="xmlfile",
                              help="name of Landsat XML file",
                              metavar="FILE")
            parser.add_option("-s", "--process_sr", type="string",
                              dest="process_sr",
                              help=("process the surface reflectance products;"
                                    " True or False (default is True) "
                                    " If False, then processing will halt"
                                    " after the TOA reflectance products are"
                                    " complete."))
            parser.add_option("--usebin",
                              dest="usebin", default=False,
                              action="store_true",
                              help=("use BIN environment variable as the"
                                    " location of LEDAPS apps"))
            parser.add_option("-l", "--logfile",
                              type="string", dest="logfile",
                              help="name of optional log file",
                              metavar="FILE")
            (options, args) = parser.parse_args()

            # validate the command-line options
            usebin = options.usebin    # should $BIN directory be used
            logfile = options.logfile  # name of the log file
            xmlfile = options.xmlfile  # name of the XML file
            if xmlfile is None:
                parser.error("missing xmlfile command-line argument")
                return ERROR
            process_sr = options.process_sr  # process SR or not
            if process_sr is None:
                process_sr = "True"  # If not provided, default to True

        # open the log file if it exists; use line buffering for the output
        log_handler = None
        if logfile is not None:
            log_handler = open(logfile, 'w', buffering=1)
        logger.info('LEDAPS processing of Landsat XML file: {0}'
                    .format(xmlfile))

        # should we expect the lnd* applications to be in the PATH or in the
        # BIN directory?
        if usebin:
            # get the BIN dir environment variable
            bin_dir = os.environ.get('BIN')
            bin_dir = bin_dir + '/'
            logger.info('BIN environment variable: {0}'.format(bin_dir)
        else:
            # don't use a path to the lnd* applications
            bin_dir = ""
            logger.info('LEDAPS executables expected to be in the PATH')

        # make sure the XML file exists
        if not os.path.isfile(xmlfile):
            logger.error('Error: XML file does not exist or is not accessible'
                         ': {0}'.format(xmlfile))
            return ERROR

        # parse the XML filename, strip off the .xml
        # use the base XML filename and not the full path.
        base_xmlfile = os.path.basename(xmlfile)
        xml = re.sub('\.xml$', '', base_xmlfile)
        logger.info('Processing XML basefile: {0}'.format(xml))

        # get the path of the XML file and change directory to that location
        # for running this script.  save the current working directory for
        # return to upon error or when processing is complete.  Note: use
        # abspath to handle the case when the filepath is just the filename
        # and doesn't really include a file path (i.e. the current working
        # directory).
        mydir = os.getcwd()
        xmldir = os.path.dirname(os.path.abspath(xmlfile))
        if not os.access(xmldir, os.W_OK):
            logger.error('Path of XML file is not writable: {0}.'
                         '  LEDAPS needs write access to the XML directory.'
                         .format(xmldir))
            return ERROR
        logger.info('Changing directories for LEDAPS processing: {0}'
                    .format(xmldir))
        os.chdir(xmldir)

        # run LEDAPS modules, checking the return status of each module.
        # exit if any errors occur.
        cmdstr = "%slndpm %s" % (bin_dir, base_xmlfile)
        logger.debug('DEBUG: lndpm command: {0}'.format(cmdstr))
        (status, output) = commands.getstatusoutput(cmdstr)
        logger.info(output)
        exit_code = status >> 8
        if exit_code != 0:
            logger.error('Error running lndpm.  Processing will terminate.')
            os.chdir(mydir)
            return ERROR

        cmdstr = "%slndcal lndcal.%s.txt" % (bin_dir, xml)
        logger.debug('DEBUG: lndcal command: {0}'.format(cmdstr))
        (status, output) = commands.getstatusoutput(cmdstr)
        logger.info(output)
        exit_code = status >> 8
        if exit_code != 0:
            logger.error('Error running lndcal.  Processing will terminate.')
            os.chdir(mydir)
            return ERROR

        if process_sr == "True":
            cmdstr = "%slndsr lndsr.%s.txt" % (bin_dir, xml)
            logger.debug('DEBUG: lndsr command: {0}'.format(cmdstr))
            (status, output) = commands.getstatusoutput(cmdstr)
            logger.info(output)
            exit_code = status >> 8
            if exit_code != 0:
                msg = 'Error running lndsr.  Processing will terminate.'
                logIt(msg, log_handler)
                os.chdir(mydir)
                return ERROR

            cmdstr = "%slndsrbm.ksh lndsr.%s.txt" % (bin_dir, xml)
            logger.debug('DEBUG: lndsrbm command: {0}'.format(cmdstr))
            (status, output) = commands.getstatusoutput(cmdstr)
            logger.info(output)
            exit_code = status >> 8
            if exit_code != 0:
                msg = 'Error running lndsrbm.  Processing will terminate.'
                logIt(msg, log_handler)
                os.chdir(mydir)
                return ERROR

        # successful completion.  return to the original directory.
        os.chdir(mydir)
        logger.info('Completion of LEDAPS.')
        if logfile is not None:
            log_handler.close()
        return SUCCESS

# ##### end of Ledaps class #####

if __name__ == "__main__":
    sys.exit(Ledaps().runLedaps())
