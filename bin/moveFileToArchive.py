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
Script to safely copy a file or files to the configure default archive
file hierarchy (if absent), confirm the MD5 sum on the transferred file,
and delete the original file from the core repository filesystem. The
repository database is updated to record this change.
'''

from osqutil.config import Config
from osqutil.setup_logs import configure_logging

from logging import INFO, WARNING, DEBUG
LOGGER = configure_logging(level=WARNING)
CONFIG = Config()

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import ArchiveLocation
from osqpipe.pipeline.archive import ArchiveManager

################################################################################
if __name__ == '__main__':

  from argparse import ArgumentParser

  # Get list of Archives.
  archives = ArchiveLocation.objects.all()
  arks     = set([ a.name for a in archives ])
  defarch  = unicode(CONFIG.default_archive)
  if defarch not in arks:
    raise StandardError("Configured default archive not found in repository.")

  # Put the default archive in first place, and sort the rest alphabetically.
  arks = [ defarch ] + sorted(list(arks - set([defarch])))

  PARSER = ArgumentParser(\
    description='A script for safe copy of repository files to the default'
    + ' file archive and deletion of the originals (after a set time) in repository.')

  PARSER.add_argument('files', metavar='<filenames>', type=str, nargs='*',
                      help='The name of the file or files to archive. The files'
                      + ' must have been recorded in the repository for the script'
                      + ' to function correctly.')

  PARSER.add_argument('-a', '--all', dest='filetype', type=str, required=False,
                      help='Archive all files of a given file type.')

  PARSER.add_argument('-A', '--archive', dest='archive_name', type=str,
                      choices=arks, default = defarch,
                      help='Name of the archive where the files should be'
                      + ' copied (default=%s).' % defarch)

  PARSER.add_argument('-c', '--copy', dest='copy_only', action='store_true',
                      help='Copy file(s) to archive without regisering them as archived.')

  PARSER.add_argument('-C', '--copy_wait_archive', dest='copy_wait_archive', action='store_true',
                      help='Force archiving in following 3 stages: 1) copy files'
                      + ' to the archive, 2) wait for 5 minutes for the file system'
                      + ' to pick up the existence of the files, 3) checks md5sums'
                      + ' and register files as  archived. The option was implemented'
                      + ' to overcome a feature of CRI Archive which occasionally'
                      + ' reports files missing even after 30s since creation.')

  PARSER.add_argument('-d', '--force_delete', dest='force_delete', action='store_true',
                      help='Force source deletion even if less than N days have'
                      + ' passed since archiving.')

  PARSER.add_argument('-m', '--force_md5_check', dest='force_md5_check', action='store_true',
                      help='Force md5 check on archived file. Adds archiving'
                      + ' information in case it is missing.')

  PARSER.add_argument('-f', '--force_overwrite', dest='force_overwrite', action='store_true',
                      help='Force overwrite for files already archived.')

  PARSER.add_argument('--days-delay', dest='archive_lag', type=int, default=7,
                      help='The number of days to wait before archiving a new file (default=7).')

  PARSER.add_argument('--debug', dest='debug', action='store_true',
                      help='Set logging level to DEBUG.')

  ARGS = PARSER.parse_args()

  if ARGS.debug:
    LOGGER.setLevel(DEBUG)

  ARCHIVER = ArchiveManager(filetype          = ARGS.filetype,
                            archive           = ARGS.archive_name,
                            copy_only         = ARGS.copy_only,
                            copy_wait_archive = ARGS.copy_wait_archive,
                            force_delete      = ARGS.force_delete,
                            force_md5_check   = ARGS.force_md5_check,
                            force_overwrite   = ARGS.force_overwrite,
                            archive_lag       = ARGS.archive_lag)

  ARCHIVER.run_archival(ARGS.files)

