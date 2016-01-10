#!/usr/bin/env python
#
# $Id$

'''Script to drop alignments from the repository, given a library
code, facility and lane number on the command line. Note that files
are unaffected.'''

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
def drop_aln(libcode, facility, lanenum):

  '''Function drops all alignments for a given flowcell lane and
  library code.'''
  
  lanenum = int(lanenum)

  try:
    lane = Lane.objects.get(library__code=libcode,
                            facility__code=facility,
                            lanenum=lanenum)
  except Lane.DoesNotExist, err:
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

###############################################################################

if __name__ == '__main__':

  (LIBCODE, FACILITY, LANENUM) = sys.argv[1:]
  drop_aln(LIBCODE, FACILITY, LANENUM)

