#!/usr/bin/env python
#
# $Id$
#
# Written by Margus
# Date: 2012.10.23
#
## Dependencies:
#
# - Password free login to target scp server
# - scp
# - sed
# - /home/fnc-odompipe/software/external/bin/xsltproc
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

from .setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging('fetch_mga')

TEST_MODE = False

def fetch_mga (flowcell, flowlane, destination, nameprefix):
  """Fetches MGA report from Genologics LIMS. Returns PDF report."""
  mgafiles = []
  flowlane = int(flowlane)

  # start logging
  if TEST_MODE:
    LOGGER.setLevel(DEBUG)
  else:
    LOGGER.setLevel(INFO)

  # install lims
  lims = Lims()
  if not lims.running():
    LOGGER.error("Remote LIMS access broken... cannot continue.")
    sys.exit("LIMS not running.")

  # search lims for a lane on flowcell
  lims_fc = lims.load_flowcell(flowcell)
  if TEST_MODE:
    lims_fc.dump()
  lane = lims_fc.get_lane(flowlane)

  # Retrieve MGA file. See upstream_lims module for supported file type
  # strings (e.g. 'LANE_MGA').
  files = lane.lims_files('LANE_MGA')
  if len(files) == 0:
    LOGGER.info("No files to retrieve for %s_%d", flowcell, lane.lane)
  for lfile in files:
    filename = lfile.uri.split('/')[-1]
    LOGGER.debug("MGA PDF File: %s", filename)
    if destination:
      local_pdf = os.path.join(destination, ("%s.pdf" % nameprefix))
    else:
      local_pdf = nameprefix + ".pdf"

    # Use requests module to download file_uri to the local_pdf
    # file. We may also want the XML file if the LIMS API supplies
    # them in future.
    local_pdf = lims.get_file_by_id(lfile.lims_id, local_pdf)
    if local_pdf:
      mgafiles.append(local_pdf)

  return mgafiles
