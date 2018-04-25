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

# Written by Margus
# Date: 2012.10.23
#
## Dependencies:
#
# - Password free login to target scp server
# - scp
# - sed
# - /home/fnc-odompipe/software/external/bin/wkhtmltopdf-amd64
#
#
# Todo/Known problems: lines 157-L161 XML file related image
# download. Image file names are not unique. A solution would be to
# write files into separate folder and use OP system tmp folder..
#

'''Module providing functions used to retrieve MGA file information
from the Genologics LIMS.'''

import sys
import os
import os.path

from .upstream_lims import Lims

from osqutil.utilities import call_subprocess, CalledProcessError
from osqutil.config import Config
from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG

CONFIG = Config()
LOGGER = configure_logging('fetch_mga')

TEST_MODE = False

def fetch_mga (flowcell, flowlane, destination, nameprefix, lims_fc=None):
  """Fetches MGA report from Genologics LIMS. Returns PDF report."""
  mgafiles = []
  flowlane = int(flowlane)

  # start logging
  if TEST_MODE:
    LOGGER.setLevel(DEBUG)
  else:
    LOGGER.setLevel(INFO)

  # install lims
  if lims_fc is None:
    lims = Lims()
    if not lims.running():
      LOGGER.error("Remote LIMS access broken... cannot continue.")
      sys.exit("LIMS not running.")

    # search lims for a lane on flowcell
    lims_fc = lims.load_flowcell(flowcell)

  if TEST_MODE:
    lims_fc.dump()
  lane = lims_fc.get_lane(flowlane)

  # Retrieve Lane MGA file. See upstream_lims module for supported
  # file type strings (e.g. 'LANE_MGA').
  files = lane.lims_files('LANE_MGA')
  if len(files) == 0:
    LOGGER.info("No files to retrieve for %s_%d", flowcell, lane.lane)
  for lfile in files:
    if lfile.uri.lower()[-4:] != 'html':
      continue
    LOGGER.debug("MGA HTML URI: %s", lfile.uri)
    if destination:
      local_pdf = os.path.join(destination, ("%s.pdf" % nameprefix))
    else:
      local_pdf = nameprefix + ".pdf"

    # Convert the HTML page direct to PDF for storage in the repository.
    cmd = [ 'wkhtmltopdf-amd64', lfile.uri, local_pdf ]
    try:
      call_subprocess(cmd, path=CONFIG.hostpath)
      mgafiles.append(local_pdf)
    except CalledProcessError, err:
      LOGGER.warning("Unable to download and/or convert MGA report to PDF.")

  return mgafiles
