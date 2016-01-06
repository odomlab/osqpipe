#!/usr/bin/env python

'''Script to fill in the blanks and create fastqc reports for all
lanes in the database which don't already have them.'''

import os

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane

from osqpipe.pipeline.laneqc import LaneFastQCReport
from osqutil.config import Config

CONFIG = Config()

# Run on our newest libraries first, since they're of most interest.
for lane in Lane.objects.all().order_by('library__extra__code_text_prefix','-library__extra__code_numeric_suffix'):
  if lane.laneqc_set.count() == 0:
    print "Generating report for lane %s..." % lane
    try:
      with LaneFastQCReport(lane=lane, path=CONFIG.hostpath) as rep:
        rep.insert_into_repository() # This step will take a while.
    except Exception, err:
      print "Error encountered; skipping."
      continue
    
