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

'''Very quick script to regenerate all LIMS lane summaryurl values to
point to wherever we're supposed to be looking these days. The current
default behaviour is to keep any summaryurls which point at actual web
pages, and discard *all* broken links regardless of whether we have a
replacement or not.'''

import httplib
from urlparse import urlparse

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane
from osqpipe.pipeline.upstream_lims import Lims

def get_http_status(url):

  '''Try and load a given URL and report the HTTP status code.'''

  urlparts = urlparse(url)

  try:
    conn = httplib.HTTPConnection(urlparts.hostname, urlparts.port)
    conn.request("HEAD", urlparts.path)
    return conn.getresponse().status
  except StandardError:
    return None

# This doesn't really need to be wrapped in a database transaction.
def regenerate_urls(lims):

  '''Re-derive summary URLs and update the repository database appropriately.'''

  for lane in Lane.objects.all():

    lims_fc   = lims.load_flowcell(lane.flowcell)

    if lims_fc:
      lims_lane = lims_fc.get_lane(lane.flowlane)

      if lims_lane:
        url = lims_lane.build_summary_url()

        if url:

          # Only replace urls if they don't exist or are not found.
          if not lane.summaryurl or get_http_status(lane.summaryurl) == 404:
            print("Updating lane %d, flowcell %s" % (lane.flowlane, lane.flowcell))
            lane.summaryurl = url
            lane.save()

        else:

          # Remove currently broken links from the database.
          if lane.summaryurl and get_http_status(lane.summaryurl) == 404:
            print("Removing broken link for lane %d, flowcell %s" % (lane.flowlane, lane.flowcell))
            lane.summaryurl = ''
            lane.save()

if __name__ == '__main__':

  lims = Lims()
  regenerate_urls(lims)
  
