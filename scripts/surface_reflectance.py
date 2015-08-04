#! /usr/bin/env python
import argparse
import os
import logging
import sys
import commands


def get_logger():
        '''Replicated from ESPA logger module'''
        # Setup the Logger with the proper configuration
        logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                    ' %(levelname)-8s'
                                    ' %(filename)s:%(lineno)d:'
                                    '%(funcName)s -- %(message)s'),
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.INFO)
        return logging.getLogger(__name__)


def parse_only_xml():
        '''Will only parse --xml XML_FILENAME from cmdline.

        Precondition:
            '--xml FILENAME' exists in command line arguments
        Postcondition:
            returns xml_filename
        '''
        # Try to parse out the XML so the exe can be determined
        parse_xml = argparse.ArgumentParser(add_help=False)
        parse_xml.add_argument('--xml', action='store',
                               dest='xml_filename', required=True,
                               help='Input XML metadata file',
                               metavar='FILE')
        (temp, extra_args) = parse_xml.parse_known_args()
        return temp.xml_filename


def is_landsat8(xml_filename):
        '''Reads xml_filename for satellitecode, checks if L8.
        Precondition:
            (1) satellite_code is the first 3 characters of input product id
            (2) xml_filename is supplied as an argument
        Postcondition:
            returns True if this satellite_code is in ['lc8','lo8']
        '''
        satellite_code = xml_filename[0:3]
        l8_prefixes = ['LC8', 'LO8']
        if satellite_code in l8_prefixes:
            return True
        else:
            return False


class EXECUTE_ERROR(Exception):
    '''Raised when command in execute_cmd returns with error'''
    def __init__(self, message, *args):
        self.message = message
        Exception.__init__(self, message, *args)


def execute_cmd(cmd_string):
        '''Execute a command line and return the terminal output

        Raises:
            EXECUTE_ERROR (Stdout/Stderr)

        Returns:
            output:The stdout and/or stderr from the executed command.
        '''
        (status, output) = commands.getstatusoutput(cmd_string)

        if status < 0:
            message = "Application terminated by signal [%s]" % cmd_string
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise EXECUTE_ERROR(message)

        if status != 0:
            message = "Application failed to execute [%s]" % cmd_string
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise EXECUTE_ERROR(message)

        if os.WEXITSTATUS(status) != 0:
            message = ("Application [%s] returned error code [%d]"
                       % (cmd_string, os.WEXITSTATUS(status)))
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise EXECUTE_ERROR(message)

        return output


def get_executable(isLandsat8):
    if(isLandsat8):
        return 'do_l8_sr.py'
    else:
        return 'do_ledaps.py'


def main():
    logger = get_logger()

    xml_filename = parse_only_xml()
    isLandsat8 = is_landsat8(xml_filename)

    cmd = [get_executable(isLandsat8)]
    cmd.extend(sys.argv[1:])  # Pass all arguments through
    cmd_string = ' '.join(cmd)
    try:
        logger.info('>>'+cmd_string)
        output = execute_cmd(cmd_string)

        if len(output) > 0:
            logger.info("\n{0}".format(output))
    except EXECUTE_ERROR:
        logger.exception(('Error running {0}.'
                          'Processing will terminate.'
                          ).format(os.path.basename(__file__)))
        raise  # Re-raise so exception message will be shown.

if __name__ == '__main__':
    main()

