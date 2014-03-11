#!/usr/bin/env python
#
# $Id$

'''Given a library code, facility and lane number from the repository,
remove all alignments and their associated files from the
repository.'''

import sys
import os
import os.path
import logging
import getopt

from osqpipe.models import Project, Library, Lane, Alignment, Filetype, Facility
from django.db import transaction

from osqpipe.pipeline.config import Config

CONFIG   = Config()
FILE_TYPE = None
DUMP_LANE = False
TEST_MODE = False

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

###############################################################################

def get_options(argv):
  '''Process command-line arguments.'''
  global FILE_TYPE, DUMP_LANE, TEST_MODE
  try:
    (opts, args) = getopt.gnu_getopt(argv[1:], "t",
                                     ["filetype=", "deletelane", "test"])
  except getopt.GetoptError:
    sys.exit("Incorrect script arguments.")
  for (key, val) in opts:
    if key == "--filetype":
      FILE_TYPE = val
    if key == "--deletelane":
      DUMP_LANE = True
    if key in ("-t", "--test"):
      TEST_MODE = True
      LOGGER.setLevel(logging.DEBUG)
  return args

def unlink_file(path):
  '''Wrapper for os.unlink which knows about the dual possibilities of
  gzipped or uncompressed files present in the repository. Note that
  we catch failure to delete and warn about it, so that e.g. we can
  move the bam file out of the repository, remove the alignment and
  reload said alignment from the bam file using new pipeline code.'''
  if os.path.exists(path):
    LOGGER.info("removing '%s'" % (path,))
    if not TEST_MODE:
      try:
        os.unlink(path)
      except OSError, err:
        LOGGER.warning("Unable to delete file %s: %s", path, err)
  else:
    comp = path + CONFIG.gzsuffix
    if os.path.exists(comp):
      LOGGER.info("removing '%s'" % (comp,))
      if not TEST_MODE:
        try:
          os.unlink(comp)
        except OSError, err:
          LOGGER.warning("Unable to delete compressed file %s: %s", comp, err)
    else:
      LOGGER.warning("missing file: '%s' (and '%s')" % (path, comp))

@transaction.commit_on_success
def delete_alignment(aln, deltype):
  '''Delete an alignment from the repository, including all associated
  files.'''
  # get alnfiles

  for alnfile in aln.alnfile_set.all():
    # remove the files and drop the record from the database.
    if deltype is None or deltype == alnfile.filetype:
      LOGGER.info("Dropping alignment file '%s' from repository",
                  alnfile.filename)
      if not TEST_MODE:
        unlink_file(alnfile.repository_file_path)
        alnfile.delete()

  # Finally, delete the actual alignment record.
  if not TEST_MODE:
    aln.delete()

@transaction.commit_on_success
def delete_lane(lane):
  '''Delete a given lane from the repository, including all associated
  files.'''

  for lfile in lane.lanefile_set.all():
    # remove the files and drop the record from the database.
    LOGGER.info("Dropping lane file '%s' from repository", lfile.filename)
    if not TEST_MODE:
      unlink_file(lfile.repository_file_path)
      lfile.delete()

  # remove the lane from the db
  LOGGER.info("Dropping lane record from the repository: ID=%s", lane.id)
  if not TEST_MODE:
    lane.delete()

def delete_laneAlignments(libcode, facility, lanenum, deltype=None):
  '''Iterate over all alignments for a given lane and delete them all.'''

  if deltype: # A file type has been specified
    deltype = Filetype.objects.get(code=deltype)
  lane = Lane.objects.get(library__code=libcode,
                          facility__code=facility,
                          lanenum=int(lanenum))
  for aln in lane.alignment_set.all():
    delete_alignment(aln, deltype)
  if DUMP_LANE:
    delete_lane(lane)

###############################################################################

if __name__ == '__main__':

  LOGGER.setLevel(logging.INFO)
  (LIBCODE, FACILITY, LANENUM) = get_options(sys.argv)

  delete_laneAlignments(LIBCODE, FACILITY, LANENUM, FILE_TYPE)
