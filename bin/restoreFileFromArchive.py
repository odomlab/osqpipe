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

'''
Script to copy a file stored in the archive back to the main repository
file tree. Once the file copy operation is complete the repository
database is updated to reflect the change. No changes are made to the
archive.
'''

from osqutil.config import Config
from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)
CONFIG = Config()

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.pipeline.archive import ArchiveManager

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
