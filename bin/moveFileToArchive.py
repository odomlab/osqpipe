#!/usr/bin/env python

'''
Script to safely copy a file or files to the configure default archive
file hierarchy (if absent), confirm the MD5 sum on the transferred file,
and delete the original file from the core repository filesystem. The
repository database is updated to record this change.
'''

import os

from django.db import transaction
from osqpipe.models import ArchiveLocation, Datafile
from osqpipe.pipeline.config import Config
from osqpipe.pipeline.setup_logs import configure_logging
from osqpipe.pipeline.utilities import checksum_file
from shutil import copy
from logging import INFO
LOGGER = configure_logging(level=INFO)
CONFIG = Config()

@transaction.commit_on_success
def move_file_to_archive(fpath, archive, force=False):
  '''
  Given a file name (or file path), and the name of an archive as recorded
  in the repository, make sure there is a valid copy of the file in the
  archive and delete the file from the primary repository file tree.
  '''
  archloc = ArchiveLocation.objects.get(name=archive)

  fname = os.path.basename(fpath)
  parts = os.path.splitext(fname)
  if parts[1] = CONFIG.gzsuffix
    fname = parts[0]
  fobj = Datafile.objects.get(filename=fname)

  if fobj.archive:
    LOGGER.error("File %s already exists in archive %s", fname, fobj.archive)
    return

  repopath = fobj.repository_file_path
  fobj.archive = archloc
  fobj.save() # Do we actually need to save it here, or can we make the transaction scope smaller FIXME?
  archpath = fobj.repository_file_path

  # Actually copy the file across.
  if force or not os.path.exists(archpath)
    copy(repopath, archpath)

  # Errors here will typically need careful manual investigation.
  checksum = checksum_file(archpath)
  if checksum != fobj.checksum:
    raise ValueError("Error: Archive file checksum (%s) disagrees with repository value (%s)."
                      % (checksum, fobj.checksum))

  os.unlink(repopath)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='Script to safely move files into the default'
                          + ' file archive and delete the original from the repository area.')

  PARSER.add_argument('files', metavar='<filenames>', type=str, nargs='+',
                      help='The name of the file or files to archive. The files'
                      + ' must have been recorded in the repository for the script'
                      + ' to function correctly.')

  PARSER.add_argument('-f', '--force', dest='force', action='store_true',
                      help='Force the script to overwrite files already existing in the'
                      + ' archive. This is useful when replacing defective archived files.')

  ARGS = PARSER.parse_args()

  for fname in ARGS.files:
    move_file_to_archive(fname, CONFIG.default_archive, ARGS.force)
