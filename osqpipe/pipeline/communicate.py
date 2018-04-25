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
import os
import django
from datetime import date

from osqutil.config import Config
from osqutil.utilities import call_subprocess, parse_repository_filename, run_in_communication_host
from osqpipe.models import Lane, Status, Library, Facility, Machine
from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

class CommunicateStatus(object):
    '''This class is used for communicating lane status back to repository. Communication may happen either directly or over ssh.'''
    
    def __init__(self):

        django.setup()
        
        self.conf = Config()
        
    def setLaneStatus(self, lane_id, status):
        try:
            status_obj = Status.objects.get(code=status)
        except Status.DoesNotExist, _err:
            raise SystemExit("Status '%s' does not exist!\n" % status)
        try:
            lane = Lane.objects.get(id=lane_id)
        except Lane.DoesNotExist, _err:
            raise SystemExit("No lane with id %d!\n" % lane_id)
        lane.status = status_obj
        lane.save()
        LOGGER.info("Lane %d status changed to \'%s\'" % (lane_id, lane.status.code))


    # NB! Following function is currently not used and needs work on Machine object side!
    def setLaneStatusByFlowcell(self, library, flowcell, flowlane, facility, status):
        '''Changes status of the lane matching the library/flowcell/flowlane combination.'''
        # 1. Get lane_id matching library/flowcell/flowlane
        try:
            lib = Library.objects.search_by_name(library)
        except Library.DoesNotExist, _err:
            raise SystemExit("No library %s in repository." % library)
        # 2. Find lane and create one if not found
        try:
            lane = Lane.objects.get(flowcell=flowcell, flowlane=flowlane, library=lib)
            # 3. Set lane status
            self.setLaneStatus(int(lane.id), status)
        except Lane.DoesNotExist, _err:
            LOGGER.info("No lane with flowcell=%s, flowlane=%s, library=%s. Creating new lane." % (flowcell, flowlane, library))
            # 3. Create new lane and set status at the same time.
            # NB! Some values set here on default e.g. rundate, info about parity from lib.paired, machine etc. are set to the best of knowledge
            #     and expected to be corrected later by other processes.
            try:
                facobj = Facility.objects.get(code=facility)
            except Facility.DoesNotExist, _err:
                raise SystemExit("No facility %s in repository." % facility)
            if facility == 'EDG':
                try:
                    machine_obj = Machine.objects.get(code='EdinburghX10')                    
                except Machine.DoesNotExist, _err:
                    raise SystemExit("No machine %s in repository." % 'EdinburghX10')
            else:
                try:
                    machine_obj = Machine.objects.get(code='Unknown')
                except Machine.DoesNotExist, _err:
                    raise SystemExit("No machine %s in repository." % 'Unknown')
            try:
                status_obj = Status.objects.get(code=status)
            except Status.DoesNotExist, _err:
                raise SystemExit("Status '%s' does not exist!\n" % status)
                
            lane = Lane(library  = lib,
                        flowcell = flowcell,
                            
                        rundate  = '2008-01-01', # Default. The value should 
                        
                        lanenum  = Lane.objects.next_lane_number(lib),
                        flowlane = int(flowlane),
                        paired = lib.paired,

                        facility = facobj,
                        seqsamplepf = '',
                        seqsamplebad = '',

                        failed = False, 
                        
                        machine = machine_obj,
                        status  = status_obj)
            lane.save()
        
    def setLaneStatusByFilename(self, fname, status):
        '''Changes status of the lane matching the library/facility/lanenum info extracted from the file name.'''
        # 1. Get lane_id matching library/facility/lanenum
        (library, facility, lanenum, pipeline) = parse_repository_filename(fname)
        # 2. Find lane and create one if not found
        try:                
            lib = Library.objects.search_by_name(library)
        except Library.DoesNotExist, _err:
            raise SystemExit("No library %s in repository." % library)
        try:
            facobj = Facility.objects.get(code=facility)
        except Facility.DoesNotExist, _err:
            raise SystemExit("No facility %s in repository." % facility)
        try:                
            lane = Lane.objects.get(library=lib, facility = facobj, lanenum = lanenum)
            # 3. Set lane status
            self.setLaneStatus(int(lane.id), status)
        except Lane.DoesNotExist, _err:
            raise SystemExit("No lane code=%s, facility=%s, lanenum=%d." % (lib.code, facility, lanenum))
            
if __name__ == '__main__':

    run_in_communication_host(sys.argv)

    import argparse
    
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
