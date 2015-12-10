#!/usr/bin/env python

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
  
