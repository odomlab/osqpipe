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

'''A script to add files to alignments which already exist in the
repository.'''

import sys
import os.path
from datetime import date
from shutil import move

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqutil.utilities import parse_repository_filename, checksum_file, set_file_permissions
from osqpipe.models import Filetype, Lane, Alignment, Alnfile, Facility, Library
from django.db import transaction
from osqutil.config import Config

CONFIG = Config()

@transaction.atomic
def _save_file_to_database(fname, aln, chksum):
  '''
  Transaction-managed smallest unit of work that we can do with the
  database to save a file to a given Alignment.
  '''
  aln = Alignment.objects.get(id=aln.id) # Reload passed object within transaction.
  filetype = Filetype.objects.guess_type(fname)
  LOGGER.debug("Found filetype: %s", filetype)
  basefn = os.path.split(fname)[1]
  LOGGER.debug("basefn: '%s'", basefn)
  fnparts = os.path.splitext(basefn)
  if fnparts[1] == CONFIG.gzsuffix:
    basefn = fnparts[0]
  LOGGER.debug("basefn: '%s'", basefn)
  afile = Alnfile.objects.create(filename=basefn,
                                 checksum=chksum, filetype=filetype,
                                 alignment=aln)

  # Move files to permanent locations.
  destname = afile.repository_file_path
  LOGGER.debug("Moving %s to %s", fname, destname)
  move(fname, destname)
  set_file_permissions(CONFIG.group, destname)

  LOGGER.info("Added '%s' to '%s'", fname, aln.lane.library.code)

def append(fname, library=None, facility=None, lanenum=None, genome=None):
  '''
  Given a filename, figure out where it belongs and load it into
  the repository. Additional hints may be provided.
  '''
  LOGGER.info("Processing %s", fname)

  argcheck = [ x is not None for x in (library, facility, lanenum) ]
  if any(argcheck):
    if not all(argcheck):
      raise ValueError("Either use filename on its own, or all of"
                       + " the following: library, facility, lanenum.")
  else:
    (library, facility, lanenum, _pipeline) = parse_repository_filename(fname)

  try:
    library = Library.objects.get(code=library)
  except Library.DoesNotExist, err:
    raise StandardError("No library found with code %s" % (library,))
  facobj = Facility.objects.get(code=facility)
  lanelist = library.lane_set.filter(facility=facobj, lanenum=lanenum)
  if lanelist.count() == 0:
    LOGGER.error("Could not find lane for '%s'", fname)
  elif lanelist.count() > 1:
    LOGGER.error("Found multiple lanes for '%s': %s",
                 fname, ", ".join([x.id for x in lanelist]))
  else:
    lane = lanelist[0]

    if genome is None:
      alns = lane.alignment_set.all()
    else:
      alns = lane.alignment_set.filter(genome__code=genome)

    if alns.count() > 1:
      LOGGER.error("Too many alignments for lane '%s'; consider supplying a genome code.", lane.id)
    else:
      aln = alns[0]

      # Calculate checksum outside the db transaction, since it takes
      # a long time.
      chksum = checksum_file(fname)
      _save_file_to_database(fname, aln, chksum)

###########################################################

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(
    description='Append a file to an alignment already existing in the'
    + ' repository. Hints may be provided as to which alignment should'
    + ' be used; otherwise the script will attempt to guess based on'
    + ' the filename.')

  PARSER.add_argument('-f', '--file', dest='file', type=str, required=True,
                      help='The name of the file to append.')
  PARSER.add_argument('-l', '--library', dest='library', type=str, required=False,
                      help='The name of the library.')
  PARSER.add_argument('-p', '--facility', dest='facility', type=str, required=False,
                      help='The facility code (e.g., CRI, SAN).')
  PARSER.add_argument('-n', '--lanenum', dest='lanenum', type=int, required=False,
                      help='The lane number.')
  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=False,
                      help='The genome used in the alignment.')

  ARGS = PARSER.parse_args()

  append(fname    = ARGS.file,
         library  = ARGS.library,
         facility = ARGS.facility,
         lanenum  = ARGS.lanenum,
         genome   = ARGS.genome)
