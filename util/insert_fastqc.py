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

'''Script to process a set of pre-existing fastqc reports into the
repository. Run this script in a directory containing library-specific
subdirectories.'''

import os

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane, Lanefile

from osqpipe.pipeline.laneqc import LaneFastQCReport
from osqutil.utilities import parse_repository_filename
from osqutil.config import Config

CONFIG = Config()

for d in os.listdir('.'):
  if os.path.isdir(d):
    for r in os.listdir(d):
      if os.path.isdir(os.path.join(d,r)):

        (lib, fac, lane, _pipeline) = parse_repository_filename(r)

        if lib is None:
          print "Cannot parse dir name: %r. Skipping."
          continue
          
        laneid = "%s_%s%02d" % (lib, fac, lane)
        try:
          lane = Lane.objects.get(library__code=lib,
                                  facility__code=fac,
                                  lanenum=lane)
        except Lane.DoesNotExist, err:
          print "Lane not found in DB: %s. Skipping." % laneid
          continue
        if lane.laneqc_set.count() > 0:
          print "Lane already has QC: %s. Skipping." % laneid
          continue

        lanefiles = Lanefile.objects.filter(lane=lane,
                                            filetype__code='fq')
        fns = [ x.repository_file_path for x in lanefiles ]

        with LaneFastQCReport(target=lane,
                              path=CONFIG.hostpath,
                              workdir=os.path.abspath(d)) as rep:
  
          print "Postprocessing files for %s..." % laneid
          rep.postprocess_results(fns)

          print "Inserting %s into database..." % laneid
          rep.insert_into_repository()
        
