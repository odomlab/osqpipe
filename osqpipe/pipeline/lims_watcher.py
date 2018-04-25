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

'''Code to test for new, previously unseen FlowCell IDs in the
SolexaLIMS back-end MySQL database. Original version shamelessly
stolen from Gord Brown's runNotifier script.'''

from osqutil.utilities import munge_cruk_emails
from .upstream_lims import Lims, get_lims_run_history
from osqutil.config import Config

from ..models import Lane, Library, Status, Facility, User

################################################################################

from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging('lims_watcher')

################################################################################

class LimsWatcher(object):

  '''Class used to poll the Genologics LIMS REST API and retrieve a
  list of Run IDs ready for downstream processing. The associated
  lanes are created in the repository (assuming that the libraries are
  available), and returned to the caller. Note that typically the
  caller will want to re-query to retrieve all ready lanes, not just
  those found in the last sweep.'''

  __slots__ = ('conf', 'missing_libraries', 'user_emails', 'lims')

  def __init__(self, lims=None, debug=False):
    if debug:
      LOGGER.setLevel(DEBUG)
    else:
      LOGGER.setLevel(INFO)
    self.conf = Config()
    self.missing_libraries = set()
    self.user_emails       = set()
    if lims is None:
      lims = Lims()
    self.lims = lims

  def find_ready_lanes(self, window=3):
    '''Get a listing of FlowCell IDs ready for processing by our pipeline.'''

    people = User.objects.all()
    emails = munge_cruk_emails([x.email.lower() for x in people])

    facility = Facility.objects.get(code=self.conf.core_facility_code)
    ready    = Status.objects.get(code=self.conf.core_ready_status,
                                  authority=facility)

    # Default is a 3-day sliding window.
    root = get_lims_run_history(self.conf.lims_rest_uri, window)

    newlanes = []

    for run_elem in root.findall('.//run'):

      # If a run has a file, it must have a runFolder.
      run_id = run_elem.find('./runFolder').text
      lims_fc = self.lims.load_run(run_id)
      for limslane in lims_fc.iter_lanes():
        if any(x in emails for x in limslane.user_emails):

          # Lane contains a sample we've submitted.
          self.user_emails = self.user_emails.union(limslane.user_emails)
          lanelibs = [ x.strip() for x in limslane.user_sample_id.split(',') ]
          for libcode in lanelibs:
            try:

              lane = Lane.objects.get(flowcell=lims_fc.fcid,
                                      flowlane=limslane.lane,
                                      facility=facility,
                                      library__code__iexact=libcode)
            except Lane.DoesNotExist, _err:
              # Lane is not in the database and is therefore of
              # interest.

              # Next, see if we have sufficient info to process it.
              try:
                lib = Library.objects.get(code__iexact=libcode)

              except Library.DoesNotExist, _err:

                # However, library is not yet available. Caller
                # is responsible for detecting and handling this.
                self.missing_libraries.add(libcode)
                LOGGER.warning("Library not found in repository: %s", libcode)
                continue

              # FIXME include run ID here somewhere.
              lane = Lane(flowcell=lims_fc.fcid,
                          flowlane=limslane.lane,
                          facility=facility,
                          library=lib,
                          lanenum=Lane.objects.next_lane_number(lib),
                          status=ready)
              lane.save()

            if lane.status is ready:
              newlanes.append(lane)

    return newlanes

################################################################################

def run_limswatcher(debug=False):
  '''
  Demo function to poll the lims for newly completed flow cells and
  dump a listing to stdout.
  '''
  watcher  = LimsWatcher(debug=debug)
  newlanes = watcher.find_ready_lanes()

  # Brief demo output.
  if len(newlanes) > 0:
    newcells = set([ x.flowcell for x in newlanes])
    for fcid in newcells:
      lanes = [ x for x in newlanes if x.flowcell == fcid ]
      lanes = sorted(lanes, key=lambda x: x.flowlane)
      print "Flowcell %s completed" % fcid
      for lane in lanes:
        print "  lane %d: %s" % (lane.flowlane, lane.library.code)

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description="Poll the LIMS for newly complete flow cells")

  PARSER.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                      default=False, help="Generate debugging output")

  ARGS = PARSER.parse_args()
  run_limswatcher(debug=ARGS.verbose)
