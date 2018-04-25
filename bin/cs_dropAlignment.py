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
    for qc in aln.alignmentqc_set.all():
      for fobj in qc.alnqcfile_set.all():
        LOGGER.info("Dropping alnqcfile '%s'", fobj.filename)
        if not TEST_MODE:
          fobj.delete()
      LOGGER.info("Dropping alignment QC %d", qc.id)
      if not TEST_MODE:
        qc.delete()

    for fobj in aln.alnfile_set.all():
      LOGGER.info("Dropping alnfile '%s'", fobj.filename)
      if not TEST_MODE:
        fobj.delete()
    LOGGER.info("Dropping alignment %d", aln.id)
    if not TEST_MODE:
      aln.delete()

###############################################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Delete all alignments from the repository for a given demultiplexed Lane.')

  PARSER.add_argument('-l', '--library', dest='libcode', type=str, required=True,
                      help='The Library for which lanes should be deleted.')

  PARSER.add_argument('-f', '--facility', dest='facility', type=str, required=True,
                      help='The facility code for the sequencing (e.g. CRI).')

  PARSER.add_argument('-n', '--lane', dest='lanenum', type=int, required=True,
                      help='The flow lane number (using our own internal numbering, *not* the flowcell lane number).')

  ARGS = PARSER.parse_args()

  drop_aln(ARGS.libcode, ARGS.facility, ARGS.lanenum)
