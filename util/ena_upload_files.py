#!/bin/env python
#
# $id$

import os
import sys
from osqpipe.pipeline.ena import enaFileUploader

if __name__ == '__main__':

    from argparse import ArgumentParser
    
    PARSER = ArgumentParser(description='A script for uploading files to ENA.')

    PARSER.add_argument('-c', '--ena_credenials_file', dest='credentials_file', type=str, required=True, default=None,
                        help='ENA credentials file is a tab delimited file where the first column is either \'username\' or \'password\' followed by values in the second.')

    PARSER.add_argument('-d', '--submission_dir', dest='submission_dir', type=str, required=False, default=None,
                        help='Submission directory (used for saving md5 files for XML uploads)')
    
    PARSER.add_argument('files', metavar='N', type=str, nargs='+',
                        help='Files with their full path.')
    
    ARGS = PARSER.parse_args()

    efu = enaFileUploader(ARGS.credentials_file, ARGS.files, ARGS.submission_dir)
