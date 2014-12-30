#!/usr/bin/env python
#
# $Id$

'''Script to fetch fastq files from the LIMS.'''

from osqpipe.pipeline.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

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
