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
from logging import INFO
LOGGER = configure_logging(level=INFO)
CONFIG = Config()

def transfer_over_scp(f1, f2, port=22, user=None, attempts = 1, sleeptime = 2):
  '''A wrapper for scp allowing multipe attempts for the transfer in case of recoverable error.'''
  attempt = 0
  known_exit_messages = ['No such file or directory','Failed to add the host to the list of known hosts','Operation not permitted']

  cmd = ''
  if user is None:
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
  '''Create folder in foreign host over ssh.'''

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
  

def get_files_for_filetype(filetype):
  '''Returns file names of all files for a given file type. Currently, only two file types, fq and bam, are supported.'''

  files = []

  if filetype == 'fq':
    files = Lanefile.objects.filter(filetype__code=filetype)
  elif filetype == 'bam':
    files = Alnfile.objects.filter(filetype__code=filetype)
  else:
    raise StandardError("'%s' files not supported. Use one of the following file types: [fq, bam]" % filetype)

  LOGGER.info("%d \'%s\' files found." % (len(files), filetype))

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
def move_file_to_archive(fpath, archive, force_overwrite=False, force_delete=False, force_md5_check=False):
  '''
  Given a file name (or file path), and the name of an archive as recorded
  in the repository, make sure there is a valid copy of the file in the
  archive and delete the file from the primary repository file tree.
  Force_overwrite forces file in archive as well as archiving record in the database to be overwritten.
  Force_delete forces file in the source to be deleted even less than nr of days (specified in archive_location table) have passed since archiving.
  Force_md5_check forces md5 sum of the archived file to be checked against the record in the repository.
  '''

  archloc = ArchiveLocation.objects.get(name=archive)

  fname = os.path.basename(fpath)
  parts = os.path.splitext(fname)
  if parts[1] == CONFIG.gzsuffix:
    fname = parts[0]
  fobj = find_file(fname)
  # Setting archive=None allows retrieval of repository location of the file. NB! Do not save this as it would modify the object in db!
  archivetmp = fobj.archive
  fobj.archive = None
  repopath = fobj.repository_file_path
  fobj.archive = archivetmp

  # check if file in archive (both on db and on disk)
  alreadyInArchive = False
  if fobj.archive:
    alreadyInArchive = True
    LOGGER.warning("File %s already in archive. Date of archiving: %s." % (fname, fobj.archive_date) )
    if force_overwrite:
      fobj.archive = archloc
      fobj.archive_date = time.strftime('%Y-%m-%d')
      fobj.save() # Do we actually need to save it here, or can we make the transaction scope smaller FIXME?
      LOGGER.info("Force overwrite. Updating archive record for %s." % fname)
  else:
    fobj.archive = archloc
    fobj.archive_date = time.strftime('%Y-%m-%d') # NB! date format not tested!
    fobj.save() # Do we actually need to save it here, or can we make the transaction scope smaller FIXME?
    LOGGER.info("Creating archive record for %s." % fname)

  archpath = fobj.repository_file_path

  # Copy file to archive. In case remote host is provided, scp via remote host. Otherwise copy.
  if force_overwrite or not os.path.exists(archpath):
    
    if archloc.host and archloc.host_port and archloc.host_path and archloc.host_user:
      # TODO: fixme? The following line is a forced solution due to fobj.repository_file_path being able to deliver only repository primary location
      # and mounted archive locations of the file; while here we need the path to the file in foreign host.
      host_archdir = os.path.join(archloc.host_path, fobj.libcode)

      LOGGER.info("Creating dir for the file in archive.")
      retcode = create_foreign_dir(archloc.host, archloc.host_port, archloc.host_user, host_archdir)
      if retcode != 0:
        raise ValueError("Error: Failed to create %s in archive."
                         % (host_archdir))  
      LOGGER.info("Copying %s to the archive." % (fname))
      retcode = transfer_over_scp(repopath, '%s:%s' % (archloc.host, host_archdir), port=archloc.host_port, user=archloc.host_user)
      if retcode != 0:
        raise ValueError("Error: Failed to copy %s to archive (%s@%s:%s)." 
                         % (repopath, archloc.host_user, archloc.host, host_archdir))
    else:
      archdir = os.path.split(archpath)[0]

      if not os.path.exists(archdir):
        LOGGER.info("Creating dir for the file in archive.")
        os.makedirs(archdir)
      # The advantage of shutil.copy2 over copy is that it tries to keep the file metadata. Equivalent to 'cp -p source destination'
        LOGGER.info("Copying file to the archive.")
      t0 = time.time()
      shutil.copy2(repopath, archpath)
      time_diff = time.time() - t0
      LOGGER.info("Copying to archive completed in %d seconds." % (time_diff))
      
  elif not alreadyInArchive:
    LOGGER.info("File %s already in archive. No need to copy." % fname)
  
  # Errors here will typically need careful manual investigation.
  if (alreadyInArchive and (force_overwrite or force_md5_check)) or not alreadyInArchive:
    LOGGER.info("Comparing md5 sum of %s in archive and in repository ..." % (fname))
    checksum = checksum_file(archpath)
    if checksum != fobj.checksum:      
      raise ValueError("Error: Archive file checksum (%s) not same as in repository (%s)."
                       % (checksum, fobj.checksum))
    LOGGER.info("Md5 sum in repository and for %s are identical." % (archpath))
  
  # compute date diff
  t_date = datetime.date.today()
  pdate = str(fobj.archive_date).split('-')
  a_date = datetime.date(int(pdate[0]),int(pdate[1]),int(pdate[2]))
  date_diff = t_date - a_date
  days_diff = int(date_diff.days)

  # if file has been in Archive for long enough or force_delete, delete the source
  if force_delete or (days_diff > archloc.host_delete_timelag):
    if not force_delete:
      LOGGER.info("More than %d days passed since archiving (%d required). (Archive date=%s\tToday=%s). Removing %s" % (days_diff, archloc.host_delete_timelag, fobj.archive_date, t_date, repopath))
    else:
      if os.path.exists(archpath):
        LOGGER.info("Executing forced deletion. Archive information: date=%s file=%s. Removing %s" % (fobj.archive_date, archpath, repopath))
      else:
        raise ValueError("Error: File %s recorded to be in archive but missing on disk!" % (archpath))
    os.unlink(repopath)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='Script to safely move files into the default'
                          + ' file archive and delete the original from the repository area.')

  PARSER.add_argument('files', metavar='<filenames>', type=str, nargs='*',
                      help='The name of the file or files to archive. The files'
                      + ' must have been recorded in the repository for the script'
                      + ' to function correctly.')

  PARSER.add_argument('-a', '--all', dest='filetype',
                      help='Archive all files of a given file type.')

  PARSER.add_argument('-d', '--force_delete', dest='force_delete', action='store_true',
                      help='Force source deletion even if less than N days have passed since archiving.')

  PARSER.add_argument('-m', '--force_md5_check', dest='force_md5_check', action='store_true',
                      help='Force md5 check on archived file. Adds archiving information in case it is missing.')

  PARSER.add_argument('-f', '--force_overwrite', dest='force_overwrite', action='store_true',
                      help='Force the script to overwrite files already existing in the'
                      + ' archive. This is useful when replacing defective archived files.')

  ARGS = PARSER.parse_args()

  fnames = ARGS.files
  if ARGS.filetype:
    fnames = get_files_for_filetype(ARGS.filetype)

  for fname in fnames:
    move_file_to_archive(fname, CONFIG.default_archive, force_overwrite=ARGS.force_overwrite, force_delete=ARGS.force_delete, force_md5_check=ARGS.force_md5_check)
    # move_file_to_archive(fname,'ark', force_overwrite=ARGS.force_overwrite, force_delete=ARGS.force_delete, force_md5_check=ARGS.force_md5_check)
