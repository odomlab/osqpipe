#!/usr/bin/env python

'''
Script to safely copy a file or files to the configure default archive
file hierarchy (if absent), confirm the MD5 sum on the transferred file,
and delete the original file from the core repository filesystem. The
repository database is updated to record this change.
'''

import time

from osqpipe.models import ArchiveLocation
from osqpipe.pipeline.archive import ArchiveManager
from osqpipe.pipeline.config import Config
from osqpipe.pipeline.setup_logs import configure_logging

from logging import INFO, WARNING
LOGGER = configure_logging(level=WARNING)
CONFIG = Config()

################################################################################
if __name__ == '__main__':

  # get list of Archives
  archives = ArchiveLocation.objects.all()
  arks = []
  archive = CONFIG.default_archive
  for a in archives:
    if a.name == archive:
      arks.append('%s (default)' % a.name)
    else:
      arks.append(a.name)


  from argparse import ArgumentParser

  PARSER = ArgumentParser(\
    description='A script for safe copy of repository files to the default'
    + ' file archive and deletion of the originals (after a set time) in repository.')

  PARSER.add_argument('files', metavar='<filenames>', type=str, nargs='*',
                      help='The name of the file or files to archive. The files'
                      + ' must have been recorded in the repository for the script'
                      + ' to function correctly.')

  PARSER.add_argument('-a', '--all', dest='filetype',
                      help='Archive all files of a given file type.')

  PARSER.add_argument('-A', '--archive', dest='archive_name', default = archive,
                      help='Name of the archive where the files should be copied.'
                      + ' List of supported archives: %s' % ', '.join(arks))

  PARSER.add_argument('-c', '--copy', dest='copy_only', action='store_true',
                      help='Copy file(s) to archive without regisering them as archived.')

  PARSER.add_argument('-C', '--copy_wait_archive', dest='copy_wait_archive', action='store_true',
                      help='Force archiving in following 3 stages: 1) copy files'
                      + ' to the archive, 2) wait for 5 minutes for the file system'
                      + ' to pick up the existance of the files, 3) checks md5sums'
                      + ' and register files as  archived. The option was implemented'
                      + ' to overcome a feature of CRI Archive which occasionally'
                      + ' reports files missing even after 30s since ceration.')

  PARSER.add_argument('-d', '--force_delete', dest='force_delete', action='store_true',
                      help='Force source deletion even if less than N days have'
                      + ' passed since archiving.')

  PARSER.add_argument('-m', '--force_md5_check', dest='force_md5_check', action='store_true',
                      help='Force md5 check on archived file. Adds archiving'
                      + ' information in case it is missing.')

  PARSER.add_argument('-f', '--force_overwrite', dest='force_overwrite', action='store_true',
                      help='Force overwrite for files already archived.')

  ARGS = PARSER.parse_args()

  fnames = ARGS.files

  ARCHIVER = ArchiveManager()
  ## FIXME from here on down belongs in the ArchiveManager class.

  # Check that we are dealing with valid archive
  if ARGS.archive_name != archive and ARGS.archive_name not in arks:
    raise ValueError("""Error: Unknown archive '%s'. Exiting! """ % (ARGS.archive_name))
  LOGGER.info("""Archive is set to '%s'.""", archive)

  # In case file type was provided, override the list of files (if
  # any) provided as arguments
  if ARGS.filetype:
    fnames = ARCHIVER.get_files_for_filetype(ARGS.filetype, not_archived=True)
    LOGGER.info("Found %d non-archived files for copying.", len(fnames))

  if len(fnames) > 0:  
    # Copy files to archive
    if ARGS.copy_wait_archive or ARGS.copy_only:
      for fname in fnames:
        LOGGER.info("Copying \'%s\' to archive.", fname)
        ARCHIVER.move_file_to_archive(str(fname), ARGS.archive_name,
                              force_overwrite = ARGS.force_overwrite,
                              force_delete    = ARGS.force_delete,
                              force_md5_check = ARGS.force_md5_check,
                              copy_only=True)

    # Wait for files copied to archive to become visible in the file system
    if ARGS.copy_wait_archive:
      LOGGER.info("Copying finished. Waiting 5 minutes for file system to register copied files.")
      time.sleep(5*60)

    # Check files to archive
    failedfns = []
    if not ARGS.copy_only:
      LOGGER.warning("Archiving %d non-archived files:", len(fnames))
      for fname in fnames:
        LOGGER.warning("""Archiving '%s'.""", fname)
        failedfn = ARCHIVER.move_file_to_archive(str(fname), ARGS.archive_name,
                                          force_overwrite = ARGS.force_overwrite,
                                          force_delete    = ARGS.force_delete,
                                          force_md5_check = ARGS.force_md5_check)
        if failedfn is not None:
          failedfns.append(failedfn)
    if failedfns:
      LOGGER.error("%d failed archiving due to md5 check sum differences:", len(failedfns))
      for failedfn in failedfns:
        LOGGER.error(failedfn)

  # Check for deletion and delete primary copies of files archived long time ago.
  if ARGS.filetype:
    fnames = []
  if ARGS.force_delete or len(fnames) == 0:
    ARCHIVER.remove_primary_files(fnames, ARGS.archive_name, ARGS.filetype,
                          force_delete=ARGS.force_delete)
