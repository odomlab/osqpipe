#!/usr/bin/env python
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

'''Script to fetch fastq files from the LIMS.'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.pipeline.fetch_fastq import FQFileFetcher

###############################################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Retrieve a FASTQ file, addressed by its'
    + ' flowCell ID and lane number, from the LIMS.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Activate test mode;'
                      + ' no filesystem changes will be made.')

  PARSER.add_argument('flowCell', metavar='<flowcell>', type=str,
                      help='The flow cell ID.')

  PARSER.add_argument('flowLane', metavar='<lane>', type=str,
                      help='The flow lane ID.')

  PARSER.add_argument('destination', metavar='<dest>', type=str,
                      help='The local directory'
                      + ' into which the file should be copied.')

  ARGS = PARSER.parse_args()

  FETCHER = FQFileFetcher(destination=ARGS.destination, test_mode=ARGS.testMode)
  FETCHER.fetch(ARGS.flowCell, ARGS.flowLane)
