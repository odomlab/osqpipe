#!/usr/bin/env python

'''
Script to safely copy a file or files to the configure default archive
file hierarchy (if absent), confirm the MD5 sum on the transferred file,
and delete the original file from the core repository filesystem. The
repository database is updated to record this change.
'''

import os
import sys
import time
import datetime
import shutil
from subprocess import Popen, PIPE

from django.db import transaction
from osqpipe.models import ArchiveLocation, Lanefile, Alnfile, QCfile, Peakfile
from osqpipe.pipeline.config import Config
from osqpipe.pipeline.setup_logs import configure_logging
from osqpipe.pipeline.utilities import checksum_file
from shutil import copy
from logging import INFO, WARNING
LOGGER = configure_logging(level=WARNING)
CONFIG = Config()

################################################################################
def transfer_over_scp(f1, f2, port=22, user=None,
                      attempts = 1, sleeptime = 2):
  '''
  A wrapper for scp allowing multipe attempts for the transfer in case
  of recoverable error.
  '''
  attempt = 0
  known_exit_messages = ['No such file or directory',
                         'Failed to add the host to the list of known hosts',
                         'Operation not permitted']

  cmd = ''
  if user is not None:
    cmd = 'scp -p -o StrictHostKeyChecking=no -P %s %s %s@%s' % (port, f1, user, f2)
  else:
    cmd = 'scp -p -o StrictHostKeyChecking=no -P %s %s %s' % (port, f1, f2)

  t0 = time.time()

  while attempt < attempts:
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    (stdout, stderr) = p.communicate()
    retcode = p.wait()
    if stdout is not None:
      sys.stdout.write(stdout)
    if stderr is not None:
      sys.stderr.write(stderr)
    if retcode != 0:
      for emsg in known_exit_messages:
        if emsg in stderr:
          LOGGER.error("%s" % (stderr))
          attempt = attempts + 1
          break
      if attempt < attempts:
        attempt += 1
        LOGGER.warning("Transfer failed with following error code: \"%s\"\nTrying again (max %d times)" % (stderr, attempts-attempt))
        time.sleep(sleeptime)
    else:
      attempt = attempts + 1
  if retcode !=0:
    raise ValueError("ERROR. Failed to transfer %s. (command=\'%s\'\n" % (f1, cmd) )

  time_diff = time.time() - t0
  LOGGER.info("Copying to archive (scp) completed in %d seconds." % (time_diff))

  return retcode

def create_foreign_dir(host, port, user, folder):
  '''
  Create folder in foreign host over ssh.
  '''
  cmd = 'ssh -p %s %s@%s \'mkdir %s\'' % (port, user, host, folder)
  
  p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
  (stdout, stderr) = p.communicate()
  retcode = p.wait()
  if retcode != 0:
    if 'File exists' in stderr:      
      LOGGER.info("Directory %s@%s:%s already exists." % (user, host, folder))
      retcode = 0
    else:
      raise ValueError("ERROR. Failed to create directory in archive (cmd=\"%s\").\nSTDOUT: %s\nSTDERR: %s\n" % (cmd, stdout, stderr) )

  return retcode
  

def get_files_for_filetype(filetype, not_archived=False):
  '''
  Returns file names of all files for a given file type. Currently,
  only two file types, fq and bam, are supported. In case of
  not_archived, returns files that have not been archived.
  '''
  files = []

  if filetype == 'fq':
    if not_archived:
      files = Lanefile.objects.filter(filetype__code=filetype, archive_id__isnull=True)
    else:
      files = Lanefile.objects.filter(filetype__code=filetype)
  elif filetype == 'bam':
    if not_archived:
      files = Alnfile.objects.filter(filetype__code=filetype, archive_id__isnull=True)
    else:
      files = Alnfile.objects.filter(filetype__code=filetype)
  else:
    raise StandardError("'%s' files not supported. Use one of the following file types: [fq, bam]" % filetype)

  return files
 
def find_file(fname):

  try:
    fobj = Lanefile.objects.get(filename=fname)
  except Lanefile.DoesNotExist:
    try:
      fobj = Alnfile.objects.get(filename=fname)
    except Alnfile.DoesNotExist:
      try:
        fobj = QCfile.objects.get(filename=fname)
      except QCfile.DoesNotExist:
        try:
          fobj = Peakfile.objects.get(filename=fname)
        except Peakfile.DoesNotExist:
          raise StandardError("Datafile %s not found in repository." % fname)

  return fobj

@transaction.commit_on_success
def move_file_to_archive(fpath, archive, force_overwrite=False, force_delete=False, force_md5_check=False, copy_only=False):
  '''
  Given a file name (or file path), and the name of an archive as recorded
  in the repository, make sure there is a valid copy of the file in the
  archive and delete the file from the primary repository file tree.

  Arguments:
  force_overwrite: forces file in archive as well as archiving record in the
                   database to be overwritten.
  force_delete: forces file in the source to be deleted even less than nr of days
                (specified in archive_location table) have passed since archiving.
  force_md5_check: forces md5 sum of the archived file to be checked against the
                   record in the repository.
  copy_only: prevents the archived file to be registered as archived. This is handy
             in case we want to move the files to archive first and register them
             as archived some other time.
  '''
  archloc = ArchiveLocation.objects.get(name=archive)

  fname = os.path.basename(fpath)
  parts = os.path.splitext(fname)
  if parts[1] == CONFIG.gzsuffix:
    fname = parts[0]
  fobj = find_file(fname)
  # Setting archive=None allows retrieval of repository location of
  # the file. NB! Do not save this as it would modify the object in
  # db!
  archivetmp = fobj.archive
  fobj.archive = None
  repopath = fobj.repository_file_path
  fobj.archive = archivetmp

  # check if file in archive (both on db and on disk)
  alreadyInArchive = False
  if fobj.archive:
    alreadyInArchive = True
    LOGGER.info("File %s already in archive. Date of archiving: %s." % (fname, fobj.archive_date) )
    if force_overwrite:
      fobj.archive = archloc
      fobj.archive_date = time.strftime('%Y-%m-%d')
      # fobj.save() # Do we actually need to save it here, or can we
      # make the transaction scope smaller FIXME?
      LOGGER.warning("Force overwrite. Updating archive record for %s." % fname)
  else:
    fobj.archive = archloc
    if not copy_only:
      fobj.archive_date = time.strftime('%Y-%m-%d') # NB! date format not tested!
      # fobj.save() # Do we actually need to save it here, or can we
      # make the transaction scope smaller FIXME?
      LOGGER.info("Creating archive record for %s." % fname)
    else:
      LOGGER.info("Copying %s to archive but not recording in repository." % fname)

  archpath = fobj.repository_file_path

  # Copy file to archive. In case remote host is provided, scp via
  # remote host. Otherwise copy.
  if force_overwrite or not os.path.exists(archpath):
    
    if archloc.host and archloc.host_port and archloc.host_path and archloc.host_user:
      # TODO: fixme? The following line is a forced solution due to
      # fobj.repository_file_path being able to deliver only
      # repository primary location and mounted archive locations of
      # the file; while here we need the path to the file in foreign
      # host.
      host_archdir = os.path.join(archloc.host_path, fobj.libcode)

      LOGGER.info("Creating dir for the file in archive.")
      retcode = create_foreign_dir(archloc.host, archloc.host_port,
                                   archloc.host_user, host_archdir)
      if retcode != 0:
        raise ValueError("Error: Failed to create %s in archive."
                         % (host_archdir))  
      LOGGER.info("Copying %s to the archive." % (fname))
      retcode = transfer_over_scp(repopath, '%s:%s' % (archloc.host, host_archdir),
                                  port=archloc.host_port, user=archloc.host_user)
      if retcode != 0:
        raise ValueError("Error: Failed to copy %s to archive (%s@%s:%s)." 
                         % (repopath, archloc.host_user, archloc.host, host_archdir))
    else:
      archdir = os.path.split(archpath)[0]

      if not os.path.exists(archdir):
        LOGGER.info("Creating dir for the file in archive.")
        os.makedirs(archdir)
        # The advantage of shutil.copy2 over copy is that it tries to
        # keep the file metadata. Equivalent to 'cp -p source
        # destination'
        LOGGER.info("Copying file to the archive.")
      t0 = time.time()
      shutil.copy2(repopath, archpath)
      time_diff = time.time() - t0
      LOGGER.info("Copying to archive completed in %d seconds." % (time_diff))
      
  elif not alreadyInArchive:
    LOGGER.info("File %s already in archive. No need to copy." % fname)
  
  if copy_only:
    return None

  # Errors here will typically need careful manual investigation.
  if (alreadyInArchive and (force_overwrite or force_md5_check)) or not alreadyInArchive:
    LOGGER.info("Comparing md5 sum of %s in archive and in repository ..." % (fname))
    checksum = checksum_file(archpath)
    if checksum != fobj.checksum:      
      # raise ValueError("Error: Archive file checksum (%s) not same as in repository (%s)."
      #                  % (checksum, fobj.checksum))
      LOGGER.error("Error: Archive file checksum (%s) not same as in repository (%s). Skipping!"
                   % (checksum, fobj.checksum))
      return fname
    else:
      fobj.save() # Do we actually need to save it here, or can we make the transaction scope smaller FIXME?
    LOGGER.info("Md5 sum in repository and for %s are identical." % (archpath))
  return None

def remove_primary_files(files, archive, filetype, force_delete=False):
  '''
  Deletes primary copies of the archived files in case forced or
  archived more than specified number of days ago
  (archloc.host_delete_timelag).  In case list of files is empty,
  identifies all files that have been in archive more than specified
  number of days for a particular file type.
  '''
  # in case empty list of files (i.e. consider deleting all 
  if len(files) == 0:
    archloc = ArchiveLocation.objects.get(name=archive)
    time_threshold = datetime.datetime.now() - datetime.timedelta(days=archloc.host_delete_timelag)
    t_date = datetime.date.today()
    if filetype == 'fq':
      files = Lanefile.objects.filter(filetype__code=filetype, archive_date__lt=time_threshold)
    if filetype == 'bam':
      files = Alnfile.objects.filter(filetype__code=filetype, archive_date__lt=time_threshold)
    LOGGER.warning("Identified %d %s files archived more than %d days ago.",
                   len(files), filetype, archloc.host_delete_timelag)
  if ARGS.force_delete:
    LOGGER.warning("Attention: FORCED DELETION of primary copies for %d archived files!",
                   len(files))

  # if file has been in Archive for long enough or force_delete, delete the source
  filesDeleted = 0
  for f in files:
    fobj = find_file(f)
    archpath = fobj.repository_file_path
    fobj.archive = None
    repopath = fobj.repository_file_path
    if not force_delete:
      LOGGER.info("More than %d days passed since archiving %s. (Archive date=%s\tToday=%s).",
                  archloc.host_delete_timelag, repopath, fobj.archive_date, t_date)
      if os.path.exists(repopath):
        LOGGER.warning("Removing %s." % repopath)
    else:
      LOGGER.warning("Executing forced deletion. Archive information: date=%s file=%s. Removing %s",
                     fobj.archive_date, archpath, repopath)

    # Before deleting the file, check that the file in archive not
    # only exists but has the same md5 sum as recorded in repository
    if os.path.exists(archpath):
      checksum = checksum_file(archpath)
      if checksum != fobj.checksum:
        # raise ValueError("Error: Archive file checksum (%s) not same
        # as in repository (%s)." % (checksum, fobj.checksum))
        LOGGER.error(\
          "Error: Archive file checksum (%s) not same as in repository (%s). Can not delete the file!",
          checksum, fobj.checksum)
      else:
        if os.path.exists(repopath):
          os.unlink(repopath)
          filesDeleted += 1
    else:
      raise ValueError("Error: File %s recorded to be in archive but missing on disk!" % (archpath))
  if filesDeleted == 0:
    LOGGER.warning("No files to delete from repository.")
  else:
    LOGGER.warning("%d files removed from repository.", filesDeleted)

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

  # Check that we are dealing with valid archive
  if ARGS.archive_name != archive and ARGS.archive_name not in arks:
    raise ValueError("Error: Unknown archive \'%s\'. Exiting! " % (ARGS.archive_name))
  LOGGER.info("Archive is set to \'%s\'." % archive)

  # In case file type was provided, override the list of files (if
  # any) provided as arguments
  if ARGS.filetype:
    fnames = get_files_for_filetype(ARGS.filetype, not_archived=True)
    LOGGER.warning("Found %d non-archived files for copying." % len(fnames))

  if len(fnames) > 0:  
    # Copy files to archive
    if ARGS.copy_wait_archive or ARGS.copy_only:
      for fname in fnames:
        LOGGER.warning("Copying \'%s\' to archive." % fname)
        move_file_to_archive(str(fname), ARGS.archive_name,
                             force_overwrite=ARGS.force_overwrite,
                             force_delete=ARGS.force_delete,
                             force_md5_check=ARGS.force_md5_check,
                             copy_only=True)

    # Wait for files copied to archive to become visible in the file system
    if ARGS.copy_wait_archive:
      LOGGER.warning("Copying finished. Waiting 5 minutes for file system to register copied files.")
      time.sleep(5*60)

    # Check files to archive
    failedfns = []
    if not ARGS.copy_only:
      LOGGER.warning("Archiving %d non-archived files:" % len(fnames))
      for fname in fnames:
        LOGGER.warning("Archiving \'%s\'." % fname)
        failedfn = move_file_to_archive(str(fname), ARGS.archive_name,
                                        force_overwrite=ARGS.force_overwrite,
                                        force_delete=ARGS.force_delete,
                                        force_md5_check=ARGS.force_md5_check)
        if failedfn is not None:
          failedfns.append(failedfn)
    if failedfns:
      LOGGER.error("%d failed archiving due to md5 check sum differences:", len(failedfns))
      for failedfn in failedfns:
        LOGGER.error("%s" % failedfn)

  # Check for deletion and delete primary copies of files archived long time ago.
  if ARGS.filetype:
    fnames = []
  if ARGS.force_delete or len(fnames) == 0:
    remove_primary_files(fnames, ARGS.archive_name, ARGS.filetype,
                         force_delete=ARGS.force_delete)
