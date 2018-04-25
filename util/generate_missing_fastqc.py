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
      with LaneFastQCReport(target=lane, path=CONFIG.hostpath) as rep:
        rep.insert_into_repository() # This step will take a while.
    except Exception, err:
      print "Error encountered; skipping."
      continue
    
