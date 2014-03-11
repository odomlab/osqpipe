#!/usr/bin/env python
#
# $Id$

# 1) look up lanes from flow cell
# 2) screen out lanes that are not interesting
# 3) report lanes which are
#    a) already in the db
#    b) not in the db, but the library is known
#    c) not in the db, and no library matches
# Also report overall status of flowcell: complete, started, etc

'''Code to check the status of a given flowcell. When run as a script,
provides a handy summary of the lanes run and how they relate to the
current contents of the repository.'''

# All script code moved into our main pipeline library namespace.
from osqpipe.pipeline.flowcell import FlowCellQuery

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

###############################################################################
          
if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Retrieve FlowCell information from the LIMS.')

  PARSER.add_argument('flowCell', metavar='<flowCellID>', type=str,
                      help='The flow cell ID (required).')

  PARSER.add_argument('flowLane', metavar='<flowLaneID>', type=int, nargs='?',
                      help='The flow lane ID (optional).')

  PARSER.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                      help='Turn on verbose output.')

  ARGS = PARSER.parse_args()

  QUERY = FlowCellQuery(flowCell=ARGS.flowCell,
                        flowLane=ARGS.flowLane,
                        verbose=ARGS.verbose)
