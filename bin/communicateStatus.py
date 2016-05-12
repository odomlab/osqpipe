#!/usr/bin/env python
#
# $Id$

import sys
import os
from communicate import CommunicateStatus

if __name__ == '__main__':
    import argparse

    PARSER = argparse.ArgumentParser(
        description='communicate Lane status to repository.')

    PARSER.add_argument('--lane_id', dest='lane_id', type=int, help='Lane ID', required = True)
    PARSER.add_argument('--status', dest='status', type=str, help='Lane status. E.g. \'qc in progress\' or \'run complete\'', required = True)

    ARGS = PARSER.parse_args()

#    lane_id = 82957 #(status_id=62172) status = 'complete' # (status_id=62175)
#    status = 'qc in progress' # (status_id=62175)
    cm = CommunicateStatus()
    cm.setLaneStatus(ARGS.lane_id, ARGS.status)
