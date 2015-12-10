#!/usr/bin/env python

'''Script to process a set of pre-existing fastqc reports into the
repository. Run this script in a directory containing library-specific
subdirectories.'''

import os

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane, Lanefile

from osqpipe.pipeline.laneqc import LaneFastQCReport
from osqpipe.pipeline.utilities import parse_repository_filename
from osqpipe.pipeline.config import Config

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

        with LaneFastQCReport(lane=lane,
                              path=CONFIG.hostpath,
                              workdir=os.path.abspath(d)) as rep:
  
          print "Postprocessing files for %s..." % laneid
          rep.postprocess_results(fns)

          print "Inserting %s into database..." % laneid
          rep.insert_into_repository()
        
