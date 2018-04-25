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

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

# All script code moved into our main pipeline library namespace.
from osqpipe.pipeline.flowcell import FlowCellQuery

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

  QUERY = FlowCellQuery(flowcell_id=ARGS.flowCell,
                        flow_lane=ARGS.flowLane,
                        verbose=ARGS.verbose)
