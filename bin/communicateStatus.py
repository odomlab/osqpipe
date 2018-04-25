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

import sys
from osqutil.utilities import run_in_communication_host

if __name__ == '__main__':

    run_in_communication_host(sys.argv)

    import argparse
    from osqpipe.pipeline.communicate import CommunicateStatus
    
    PARSER = argparse.ArgumentParser(
        description='communicate Lane status to repository. Creates new lane in case lane with flowcell & flowlane & library & status is missing.')
    
    PARSER.add_argument('--status', dest='status', type=str, help='Lane status. E.g. \'qc in progress\' or \'run complete\'', required = True)
    PARSER.add_argument('--flowlane', dest='flowlane', type=str, help='flowlane of the lane', required = False)
    PARSER.add_argument('--library', dest='library', type=str, help='library (i.e. donumber)', required = False)
    PARSER.add_argument('--facility', dest='facility', type=str, help='Facility code (e.g. CRI)', required = False)
    GROUP = PARSER.add_mutually_exclusive_group()
    GROUP.add_argument('--laneid', dest='laneid', type=int, help='Lane ID')
    GROUP.add_argument('--flowcell', dest='flowcell', type=str, help='Flowcell of the lane')
    GROUP.add_argument('--filename', dest='fname', type=str, help='Filename')

    ARGS = PARSER.parse_args()

    #   lane_id = 82957 #(status_id=62172) status = 'complete' # (status_id=62175)
    #   status = 'qc in progress' # (status_id=62175)
    cm = CommunicateStatus()

    if ARGS.laneid is not None:
        cm.setLaneStatus(ARGS.laneid, ARGS.status)
    elif ARGS.flowcell:
        if ARGS.flowlane is not None and ARGS.library is not None and ARGS.facility is not None:
            cm.setLaneStatusByFlowcell(ARGS.library, ARGS.flowcell, ARGS.flowlane, ARGS.facility, ARGS.status)
        else:
            raise SystemExit("To set lane status by flowcell, flowlane, facility and library need to be specified.")
    elif ARGS.fname is not None:
        cm.setLaneStatusByFilename(ARGS.fname, ARGS.status)
