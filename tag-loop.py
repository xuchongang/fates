#!/usr/bin/env python
"""FIXME: A nice python program to do something useful.

Author: Ben Andre <andre@ucar.edu>

"""

from __future__ import print_function

import sys

if sys.hexversion < 0x02070000:
    print(70 * "*")
    print("ERROR: {0} requires python >= 2.7.x. ".format(sys.argv[0]))
    print("It appears that you are running python {0}".format(
        ".".join(str(x) for x in sys.version_info[0:3])))
    print(70 * "*")
    sys.exit(1)

#
# built-in modules
#
import argparse
import os
import subprocess
import traceback

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser

#
# installed dependencies
#

#
# other modules in this package
#

# -------------------------------------------------------------------------------
#
# User input
#
# -------------------------------------------------------------------------------

def commandline_options():
    """Process the command line arguments.

    """
    parser = argparse.ArgumentParser(
        description='FIXME: python program template.')

    parser.add_argument('--backtrace', action='store_true',
                        help='show exception backtraces as extra debugging '
                        'output')

    parser.add_argument('--debug', action='store_true',
                        help='extra debugging output')

    parser.add_argument('--config', nargs=1, required=True,
                        help='path to config file')

    options = parser.parse_args()
    return options


# -------------------------------------------------------------------------------
#
# work functions
#
# -------------------------------------------------------------------------------
def read_config_file(filename):
    """Read the configuration file and process

    """
    print("Reading configuration file : {0}".format(filename))

    cfg_file = os.path.abspath(filename)
    if not os.path.isfile(cfg_file):
        raise RuntimeError("Could not find config file: {0}".format(cfg_file))

    config = config_parser()
    config.read(cfg_file)

    return config


def write_config_file(config, filename):
    """Read the configuration file and process

    """
    print("Writing configuration file : {0}".format(filename))

    cfg_file = os.path.abspath(filename)
    if not os.path.isfile(cfg_file):
        raise RuntimeError("Could not find config file: {0}".format(cfg_file))

    with open(cfg_file, 'w') as configfile:
        config.write(configfile)


# -------------------------------------------------------------------------------
#
# main
#
# -------------------------------------------------------------------------------

def main(options):
    tags = ['clm4_5_2_r122',
            'clm4_5_2_r123',
            'clm4_5_2_r124',
            'clm4_5_2_r125',
            'clm4_5_2_r126',
            'clm4_5_2_r127',
            'clm4_5_2_r128',
            'clm4_5_3_r129',
            'clm4_5_3_r130',
            'clm4_5_3_r131',
            'clm4_5_3_r132',
            'clm4_5_3_r133',
            'clm4_5_3_r134',
            'clm4_5_3_r135',
            'clm4_5_3_r136',
            'clm4_5_3_r137',
            'clm4_5_3_r138',
            'clm4_5_3_r139',
            'clm4_5_3_r140',
            'clm4_5_3_r141',
            'clm4_5_3_r142',
            'clm4_5_3_r143',
            'clm4_5_3_r144',
            'clm4_5_3_r145',
            'clm4_5_3_r146',
            'clm4_5_3_r147',
            'clm4_5_3_r148',
            'clm4_5_3_r149',
            'clm4_5_4_r150',
            'clm4_5_4_r151',
            'clm4_5_5_r152',
            'clm4_5_6_r153',
            'clm4_5_6_r154',
            'clm4_5_6_r155',
            'clm4_5_6_r156',
            'clm4_5_6_r157',
            'clm4_5_6_r158',
            ]

    config_filename = options.config[0]
    for t in tags:
        config = read_config_file(config_filename)

        tag = "clm2/trunk_tags/{0}".format(t)
        config.set('cesm', 'tag', tag)

        write_config_file(config, config_filename)

        cmd = ['./ed-clm-git/cesm2git.py',
               '--config',
               config_filename,
               '--feelin-lucky',
               ]
        print(tag)
        subprocess.check_output(cmd, shell=False, stderr=subprocess.STDOUT)

    return 0


if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as error:
        print(str(error))
        if options.backtrace:
            traceback.print_exc()
        sys.exit(1)
