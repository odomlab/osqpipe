#!/bin/env python
#
# Copyright 2018 Odom Lab, CRUK-CI, University of Cambridge
#
# This file is part of the osqpipe python package.
#
# The osqpipe python package is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# The osqpipe python package is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the osqpipe python package.  If not, see
# <http://www.gnu.org/licenses/>.

import os
import sys
from subprocess import Popen, PIPE
import tempfile
#from logging import INFO
#LOGGER = configure_logging(level=INFO)

# By lukk01.

def read_credentials(credentials_file):
    '''Parse ENA credentials file.'''

    credentials = {'username':None, 'password':None}
    with open(credentials_file,'rb') as fh:
        for line in fh:
            line = line.rstrip('\n')
            cols = line.split('\t')
            if cols[0] in credentials:
                credentials[cols[0]] = cols[1]
        fh.close()
    for key in credentials:
        if credentials[key] is None:
            sys.stderr.write("Ill formated credentials file. No value for \'%s\'!\n\n" % key)
            sys.exit(1)

    return credentials

def create_md5_file(fn):

    md5fn = fn + '.md5'
    if os.path.isfile(md5fn):
        return md5fn

    base = os.path.basename(fn)
    temp_dir = tempfile.gettempdir()
    md5fn = os.path.join(temp_dir, base + '.md5')

    cmd = "md5sum %s > %s" % (fn, md5fn)

    sys.stdout.write("Creating %s ...\n" % md5fn)
    subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    (stdout, stderr) = subproc.communicate()
    retcode = subproc.wait()
    #if stdout is not None:
    #  sys.stdout.write(stdout)
    #if stderr is not None:
    #  sys.stderr.write(stderr)
    if retcode != 0:
        LOGGER.error("%s\nFailed to create %s.\n\n" % md5fn)
        sys.exit(1)

    return md5fn

def ena_upload_file(credentials, fn):
    '''Uses ascp to upload a single data file to ENA webin using credentials provided.'''
    
    # Note that we are uploading gzipped file.
    cmd = "ascp -QT -l300M -L- " + fn + " " + credentials['username'] + "@webin.ebi.ac.uk:."
    sys.stdout.write("Uploading %s ...\n" % fn)
    print cmd
    subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    (stdout, stderr) = subproc.communicate()
    retcode = subproc.wait()
    #if stdout is not None:
    #  sys.stdout.write(stdout)
    #if stderr is not None:
    #  sys.stderr.write(stderr)
    if retcode != 0:
        sys.stderr.write("%s\nFile transfer failed!\n\n" % stderr)
        sys.exit(1)    

def ena_upload_files(credentials_file, files):
    '''Computes md5 sums and uploads files with sums to ENA webin using credentials from credentials_file.'''

    credentials = read_credentials(credentials_file)
    
    os.environ['ASPERA_SCP_PASS'] = credentials['password']    

    for f in files:
        md5f = create_md5_file(f)
        ena_upload_file(credentials, md5f)
        os.unlink(md5f)
        ena_upload_file(credentials, f)

    sys.exit(0)
        
if __name__ == '__main__':

    from argparse import ArgumentParser
    
    # NB! The script is deliberately design to be light weight with no dependencies to other packages.
    #     The only two requirements are ascp and md5sum being present in path.
    
    PARSER = ArgumentParser(description='A script for uploading files to ENA.')

    PARSER.add_argument('-c', '--ena_credenials_file', dest='credentials_file', type=str, required=True, default=None,
                        help='ENA credentials file is a tab delimited file where the first column is either \'username\' or \'password\' followed by values in the second.')

    PARSER.add_argument('files', metavar='N', type=str, nargs='+',
                        help='Files with their full path.')
    
    ARGS = PARSER.parse_args()
    ena_upload_files(ARGS.credentials_file, ARGS.files)
