#!/usr/bin/env python
#
# $Id$

import sys
import os
import socket
from osqutil.config import Config
from osqutil.utilities import call_subprocess
from osqpipe.models import Lane, Status
from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

def get_my_ip():
    return socket.gethostbyname(socket.getfqdn())

def get_my_hostname():
    return socket.gethostname()

class CommunicateStatus(object):
    '''This class is used for communicating lane status back to repository. Communication may happen either directly or over ssh.'''
    
    def __init__(self):

        self.conf = Config()
        
        self.ssh = False
        my_name = get_my_hostname()

        # if my_host name is different from the hostname where the communication should take place, use different host name
        if  my_name != self.conf.communicationhost:
            self.ssh = True
            self.user = self.conf.user
            self.host = self.conf.communicationhost
            self.port = self.conf.communicationport

    def setLaneStatus(self, lane_id, status):
        if self.ssh:
            self.setLaneStatusOverSsh(lane_id, status)
        else:            
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
            LOGGER.info("Lane %d status changed to \'%s\'\n" % (lane_id, status))

    def getLaneStatusOverSshCommand(self, lane_id, status):

        cmd = 'communicateStatus.py --lane %d --status %s' % (lane_id, status)
        ssh_cmd = 'ssh %s@%s %s' % (self.user, self.host, cmd)

        return ssh_cmd
        
    def setLaneStatusOverSsh(self, lane_id, status):

        ssh_cmd = getLaneStatusOverSshCommand(self, lane_id, status)
        call_subprocess(ssh_cmd, path=self.conf.hostpath, shell=True)

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
