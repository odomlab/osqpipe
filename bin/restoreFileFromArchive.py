#!/usr/bin/env python

'''
Script to copy a file stored in the archive back to the main repository
file tree. Once the file copy operation is complete the repository
database is updated to reflect the change. No changes are made to the
archive.
'''

from django.db import transaction
from osqpipe.models import Datafile
from osqpipe.pipeline.config import Config
from osqpipe.pipeline.setup_logs import configure_logging
from osqpipe.pipeline.utilities import checksum_file
from shutil import copy
from logging import INFO
LOGGER = configure_logging(level=INFO)
CONFIG = Config()

@transaction.commit_on_success
def restore_file_from_archive(fpath):

  fname = os.path.basename(fpath)
  parts = os.path.splitext(fname)
  if parts[1] = CONFIG.gzsuffix
    fname = parts[0]
  fobj = Datafile.objects.get(filename=fname)

  if not fobj.archive:
    LOGGER.error("File %s is not currently registered to any archive location.", fname)

  archpath = fobj.repository_file_path
  fobj.archive = None
  fobj.save()  # FIXME is this necessary? see moveFileToArchive
  repopath = fobj.repository_file_path

  # This really shouldn't happen, and indicates manual intervention may be necessary.
  if os.path.exists(repopath):
    raise StandardError("File %s already present in repository tree." % fname)

  copy(archpath, repopath)

  checksum = checksum_file(repopath)

  # Another manual investigation type error.
  if checksum != fobj.checksum:
    raise ValueError("Restored file checksum (%s) does not agree with repository value (%s)."
                      % (checksum, fobj.checksum))

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='Script to restore files from'
                          + ' file archives back into the repository area. The'
                          + ' original archive copy remains unchanged.')

  PARSER.add_argument('files', metavar='<filenames>', type=str, nargs='+',
                      help='The name of the file or files to restore. The files'
                      + ' must have been recorded in the repository for the script'
                      + ' to function correctly.')

  ARGS = PARSER.parse_args()

  for fname in ARGS.files:
    restore_file_from_archive(fname)
