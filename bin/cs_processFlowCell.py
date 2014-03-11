#!/usr/bin/env python
#
# $Id$

'''Script to query a given flowcell ID for lanes of interest, download
the sequencing data (i.e. fastq file) and demultiplex if necessary.'''

# All script code moved into our main pipeline library namespace.
from osqpipe.pipeline.flowcell import FlowCellProcess

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

###############################################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Process a FlowCell from the LIMS.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Turn on test mode (no filesystem changes).')

  PARSER.add_argument('-n', '--nocheck', dest='checkForLibInDB',
                      action='store_false',
                      help='Deactivate the check for library in repository DB.')

  PARSER.add_argument('flowCell', metavar='<flowCellID>', type=str,
                      help='The flow cell ID (required).')

  PARSER.add_argument('flowLane', metavar='<flowLaneID>', type=int, nargs='?',
                      help='The flow lane ID (optional).')

  PARSER.add_argument('-f', '--force-primary', dest='forcePrimary',
                      action='store_true',
                      help='Start processing a flow cell which is only PRIMARY COMPLETE,'
                      + ' not SECONDARY COMPLETE. It is the user\'s responsibility to ensure'
                      + ' that the LIMS files are ready to be copied across.')

  ARGS = PARSER.parse_args()

  PROC = FlowCellProcess(testMode=ARGS.testMode,
                         checkForLibInDB=ARGS.checkForLibInDB,
                         forcePrimary=ARGS.forcePrimary)

  PROC.run(ARGS.flowCell, ARGS.flowLane)

