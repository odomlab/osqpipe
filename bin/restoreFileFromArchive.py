#!/usr/bin/env python

'''
Script to copy a file stored in the archive back to the main repository
file tree. Once the file copy operation is complete the repository
database is updated to reflect the change. No changes are made to the
archive.
'''

from osqpipe.pipeline.archive import ArchiveManager
from osqpipe.pipeline.config import Config
from osqpipe.pipeline.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)
CONFIG = Config()

################################################################################
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

  ARCHIVER = ArchiveManager()

  for fname in ARGS.files:
    ARCHIVER.restore_file_from_archive(fname)
