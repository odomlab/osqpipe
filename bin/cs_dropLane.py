#!/usr/bin/env python
#
# $Id$

'''Script to delete a lane from the repository given library code,
facility and lane number. Note that files are unaffected.'''

import sys

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane, Library, Alignment

TEST_MODE = False

###############################################################################

def drop_lane(libcode, facility, lanenum):
  '''Given a library code, facility and lane number delete all records
  relating to that lane from the repository.'''
  lanenum = int(lanenum)

  lane = Lane.objects.get(library__code=libcode,
                          facility__code=facility,
                          lanenum=lanenum)
  if lane == None:
    errtuple = (libcode, facility, lanenum)
    LOGGER.warning("Failed to load lane '%s_%s%02d'.  Quitting.", *errtuple)
    sys.exit("Lane '%s_%s%02d' not found in repository" % errtuple)

  for aln in lane.alignment_set.all():
    for fobj in aln.alnfile_set.all():
      LOGGER.info("Dropping alnfile '%s'", fobj.filename)
      if not TEST_MODE:
        fobj.delete()
    LOGGER.info("Dropping alignment %d", aln.id)
    if not TEST_MODE:
      aln.delete()

  for qc in lane.laneqc_set.all():
    for fobj in qc.qcfile_set.all():
      LOGGER.info("Dropping qcfile '%s'", fobj.filename)
      if not TEST_MODE:
        fobj.delete()
    LOGGER.info("Dropping lane QC %d", qc.id)
    if not TEST_MODE:
      qc.delete()
      
  for fobj in lane.lanefile_set.all():
    LOGGER.info("Dropping lanefile '%s'", fobj.filename)
    if not TEST_MODE:
      fobj.delete()
  LOGGER.info("Dropping lane %s %s%02d", libcode, facility, lanenum)
  if not TEST_MODE:
    lane.delete()

###############################################################################

if __name__ == '__main__':
  (LIBCODE, FACILITY, LANENUM) = sys.argv[1:]
  drop_lane(LIBCODE, FACILITY, LANENUM)
  
