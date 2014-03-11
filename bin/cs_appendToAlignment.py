#!/usr/bin/env python
#
# $Id$

'''A script to add files to alignments which already exist in the
repository.'''

import sys
import logging
import os.path
from datetime import date
from shutil import move

from osqpipe.pipeline.utilities import parse_repository_filename, checksum_file
from osqpipe.models import Filetype, Lane, Alignment, Alnfile, Facility, Library
from django.db import transaction
from osqpipe.pipeline.config import Config

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

CONFIG = Config()

@transaction.commit_on_success
def _save_file_to_database(fname, aln, chksum):
  '''
  Transaction-managed smallest unit of work that we can do with the
  database to save a file to a given Alignment.
  '''
  filetype = Filetype.objects.guess_type(fname)
  basefn = os.path.split(fname)[1]
  LOGGER.debug("basefn: '%s'" % (basefn))
  fnparts = os.path.splitext(basefn)
  if fnparts[1] == CONFIG.gzsuffix:
    basefn = fnparts[0]
  LOGGER.debug("basefn: '%s'" % (basefn))
  afile = Alnfile.objects.create(filename=basefn,
                                 checksum=chksum, filetype=filetype,
                                 description='', alignment=aln)

  # Move files to permanent locations.
  destname = afile.repository_file_path
  LOGGER.debug("mv %s %s" % (fname, destname))
  move(fname, destname)

  LOGGER.info("Added '%s' to '%s'" % (fname, aln.lane.library.code))

def append(fname):
  '''Given a filename, figure out where it belongs and load it into
  the repository.'''
  LOGGER.info(fname)
  (code, facility, lanenum, _pipeline) = parse_repository_filename(fname)
  try:
    library = Library.objects.get(code=code)
  except Library.DoesNotExist, err:
    raise StandardError("No library found with code %s" % (code,))
  facobj = Facility.objects.get(code=facility)
  lanelist = library.lane_set.filter(facility=facobj, lanenum=lanenum)
  if lanelist.count() == 0:
    LOGGER.error("Could not find lane for '%s'", fname)
  elif lanelist.count() > 1:
    LOGGER.error("Found multiple lanes for '%s': %s",
                 fname, ", ".join([x.id for x in lanelist]))
  else:
    lane = lanelist[0]
    alns = lane.alignment_set.all()
    if alns.count() > 1:
      LOGGER.error("Too many alignments for lane '%s'" % (lane.id))
    else:
      aln = alns[0]

      # Calculate checksum outside the db transaction, since it takes
      # a long time.
      chksum = checksum_file(fname)
      _save_file_to_database(fname, aln, chksum)

###########################################################

if __name__ == '__main__':

  LOGGER.setLevel(logging.INFO)
  append(sys.argv[1])
