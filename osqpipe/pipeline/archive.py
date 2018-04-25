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
Code used to manage the transfer and restoration of files to and from
various ArchiveLocations.
'''

import os
import sys
import time
import datetime
import re

from subprocess import Popen, PIPE
from shutil import copy2

from django.db import transaction
from ..models import ArchiveLocation, Lanefile, Alnfile, \
    QCfile, AlnQCfile, Peakfile, MergedAlnfile, Datafile
from osqutil.utilities import checksum_file, bash_quote

from osqutil.config import Config
from osqutil.setup_logs import configure_logging

LOGGER = configure_logging('archive')
CONFIG = Config()

################################################################################
def _archive_file_via_scp(fobj, attempts = 1, sleeptime = 2):
  '''
  A wrapper for scp allowing multiple attempts for the transfer in case
  of recoverable error.
  '''
  unrecoverable = [ 'No such file or directory',
                    'Failed to add the host to the list of known hosts',
                    'Operation not permitted' ]

  arch = fobj.archive
  if arch is None:
    raise ValueError("Attempting to transfer file to null archive location.")

  # NOTE: We may still need to double-quote spaces the destination
  # passed to scp. Double-quoting brackets ([]) does not work, though.
  host_archdir = os.path.join(arch.host_path, fobj.libcode)
  dest = os.path.join(host_archdir,
                      os.path.basename(fobj.repository_file_path))

  cmd = 'scp -p -o StrictHostKeyChecking=no'
  if arch.host_port is not None:
    cmd += ' -P %s' % str(arch.host_port)

  # Assume we're copying from the main repository to the archive.
  # Note that we need quoting of e.g. file paths containing
  # spaces.
  cmd += ' %s' % bash_quote(fobj.original_repository_file_path)

  # Double-quote the destination, as it has to get past (a) our local
  # bash, and (b) the bash on the destination machine.
  if arch.host_user is not None:
    cmd += ' %s@%s:%s' % (arch.host_user, arch.host, bash_quote(bash_quote(dest)))
  else:
    cmd += ' %s:%s' % (arch.host, bash_quote(bash_quote(dest)))

  start_time = time.time()

  while attempts > 0:
    subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    (stdout, stderr) = subproc.communicate()
    retcode = subproc.wait()
    if stdout is not None:
      sys.stdout.write(stdout)
    if stderr is not None:
      sys.stderr.write(stderr)
    if retcode != 0:
      for mesg in unrecoverable:
        if mesg in stderr:
          LOGGER.error(stderr)
          attempts = 0
          break
      attempts -= 1
      if attempts <= 0:
        break
      LOGGER.warning(\
        'Transfer failed with error code: %s\nTrying again (max %d times)',
        stderr, attempts)
      time.sleep(sleeptime)
    else:
      break

  if retcode !=0:
    raise StandardError("ERROR. Failed to transfer file. Command was:\n   %s\n"
                        % (" ".join(cmd),) )

  time_diff = time.time() - start_time
  LOGGER.info("Copying to archive (scp) completed in %d seconds.", time_diff)

def _create_archive_dir_on_host(fobj):
  '''
  Create folder in foreign host over ssh.
  '''
  arch = fobj.archive
  if arch is None:
    raise ValueError("Attempting to transfer file to null archive location.")
  folder = bash_quote(os.path.join(arch.host_path, fobj.libcode))

  # The ssh command expects double quoting (one for our bash prompt,
  # one for the server).
  cmd = [ 'ssh' ]
  if arch.host_port is not None:
    cmd += ['-p', arch.host_port ]
  if arch.host_user is not None:
    cmd += [ '%s@%s' % (arch.host_user, arch.host) ]
  else:
    cmd += [ arch.host ]
  cmd += [ 'mkdir', bash_quote(folder), '&& chmod 750', bash_quote(folder)]

  subproc = Popen(cmd, stdout=PIPE, stderr=PIPE)
  (stdout, stderr) = subproc.communicate()
  retcode = subproc.wait()
  if retcode != 0:
    if 'File exists' in stderr:
      LOGGER.info("Directory %s@%s:%s already exists.",
                  arch.host_user, arch.host, folder)
    else:
      raise StandardError(\
        "ERROR. Failed to create directory in archive"
        + (" (cmd=\"%s\").\nSTDOUT: %s\nSTDERR: %s\n"
        % (" ".join(cmd), stdout, stderr)) )

def _copy_file_to_archive(fobj):
  '''
  Copies the file to its currently configured archive location.
  '''
  archpath = fobj.repository_file_path

  if all([ getattr(fobj.archive, key) is not None
           for key in ('host', 'host_path') ]):

    # Remote archive, file transfer via scp.
    LOGGER.info("Creating dir for the file in archive.")

    # Raises StandardError upon failure to create dir.
    _create_archive_dir_on_host(fobj)

    LOGGER.info("Copying %s to the archive.", fobj)

    # Raises StandardError upon failure to transfer.
    _archive_file_via_scp(fobj)

  else:

    # Local archive, mounted on current server.
    archdir = os.path.split(archpath)[0]

    if not os.path.exists(archdir):
      LOGGER.info("Creating dir for the file in archive.")
      os.makedirs(archdir)

    # The advantage of shutil.copy2 over copy is that it tries to
    # keep the file metadata. Equivalent to 'cp -p source
    # destination'
    LOGGER.info("Copying file to the archive.")
    start_time = time.time()
    copy2(fobj.original_repository_file_path, archpath)
    time_diff = time.time() - start_time
    LOGGER.info("Copying to archive completed in %d seconds.", time_diff)

def _get_files_for_filetype(filetype, not_archived=False):
  '''
  Returns file names of all files for a given file type. Currently,
  only two file types, fq and bam, are supported. In case of
  not_archived, returns files that have not been archived.
  '''
  files = []

  if filetype == 'fq':
    if not_archived:
      files = Lanefile.objects.filter(filetype__code=filetype,
                                      archive_id__isnull=True)
    else:
      files = Lanefile.objects.filter(filetype__code=filetype)

  elif filetype == 'bam':
    if not_archived:
      files = Alnfile.objects.filter(filetype__code=filetype,
                                     archive_id__isnull=True)
    else:
      files = Alnfile.objects.filter(filetype__code=filetype)

  else:
    raise StandardError(\
      "'%s' files not recognised. Supported file types are fq or bam."
      % filetype)

  return files

def _find_file(fname):
  '''
  Simply iterate over all the different subclasses of Datafile in the
  repository, returning the appropriate file object for a given
  filename.
  '''
  fname = os.path.basename(fname)
  parts = os.path.splitext(fname)
  if parts[1] == CONFIG.gzsuffix:
    fname = parts[0]

  # One feels there should be a better way to do this using some
  # Django magic.
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
          fobj = AlnQCfile.objects.get(filename=fname)
        except AlnQCfile.DoesNotExist:
          try:
            fobj = Peakfile.objects.get(filename=fname)
          except Peakfile.DoesNotExist:
            try:
              fobj = MergedAlnfile.objects.get(filename=fname)
            except MergedAlnfile.DoesNotExist:
              raise StandardError(\
                "Datafile %s not found in repository." % fname)

  return fobj

################################################################################
class ArchiveError(StandardError):
  pass

################################################################################
class ArchiveManager(object):
  '''
  Attributes:
    archive: The name of the ArchiveLocation to be used.
    filetype: (Optional) type of files to be archived according to our
              standard rules.
    force_overwrite: Forces file in archive as well as archiving record in the
                     database to be overwritten.
    force_delete: Forces file in the source to be deleted even when less than
                  the required number of days (specified in archive_location
                  table) have passed since archiving.
    force_md5_check: Forces the MD5 sum of the archived file to be checked
                     against the corresponding record in the repository.
    copy_only: Prevents the archived file from being registered as archived.
               This is useful in cases where we want to move the files to the
               archive first, and register them as archived at a later date.
    archive_lag: A time in days which files will stay in the main repository
                 before being archived. Default behaviour is to archive as
                 soon as possible, but this is not always desirable.
  '''

  __slots__ = ('filetype', 'archive', 'copy_only', 'copy_wait_archive',
               'force_delete', 'force_md5_check', 'force_overwrite', 'archive_lag')

  def __init__(self, filetype=None, archive=CONFIG.default_archive,
               copy_only=False, copy_wait_archive=True,
               force_delete=False, force_md5_check=False,
               force_overwrite=False, archive_lag=None):
    self.filetype          = filetype
    self.archive           = ArchiveLocation.objects.get(name=archive)
    self.copy_only         = copy_only
    self.copy_wait_archive = copy_wait_archive
    self.force_delete      = force_delete
    self.force_md5_check   = force_md5_check
    self.force_overwrite   = force_overwrite
    self.archive_lag       = archive_lag

  def _set_archive_location(self, fobj, warn=True):
    '''
    Set the fobj archive location and archive_date, checking first
    that this has not already been set. We allow the logging to be
    switched off because this method is also used in a dummy run to
    set archive location temporarily. This is so that we can copy
    files without necessarily changing the database.
    '''
    already_in_archive = False
    if fobj.archive:
      already_in_archive = True
      LOGGER.debug("File %s already in archive. Date of archiving: %s.",
                   fobj, fobj.archive_date)
      if self.force_overwrite:
        fobj.archive = self.archive
        fobj.archive_date = time.strftime('%Y-%m-%d')
        if warn:
          LOGGER.warning("Force overwrite. Updating archive record for %s.",
                         fobj)
    else:
      fobj.archive = self.archive
      if not self.copy_only:
        fobj.archive_date = time.strftime('%Y-%m-%d')
        LOGGER.debug("Creating archive record for %s.", fobj)
      elif warn:
        LOGGER.warning("Copying %s to archive but not recording in repository.",
                       fobj)

    # We deliberately do not save the changes yet; that is up to the
    # calling code, which may choose to discard the archive information.
    return (fobj, already_in_archive)

  def _copy_file_to_archive_disk(self, fobj):
    '''
    Just copy the file to its designated archive disk location; make
    no changes to the database.
    '''
    # Check if file in archive (both on db and on disk).
    (fobj, _previously_archived) = self._set_archive_location(fobj, warn=False)
    archpath = fobj.repository_file_path

    # Copy file to archive. In case remote host is provided, scp via
    # remote host. Otherwise copy.
    if self.force_overwrite or not os.path.exists(archpath):
      LOGGER.warning("Copying file %s to archive location %s.",
                     fobj, fobj.archive)
      _copy_file_to_archive(fobj)

    # N.B. we do *not* want to make fobj available to the caller as we
    # are likely not within a transaction.
    return

  @transaction.atomic
  def _register_file_in_archive(self, fobj):
    '''
    Given a file name (or file path), and the name of an archive as recorded
    in the repository, make sure there is a valid copy of the file in the
    archive and delete the file from the primary repository file tree.
    '''
    LOGGER.warning("""Registering '%s' in archive.""", fobj)

    # Reload the Datafile object within this transaction.
    fobj = type(fobj).objects.get(filename=fobj.filename)
    (fobj, previously_archived) = self._set_archive_location(fobj)
    archpath = fobj.repository_file_path

    # Errors here will typically need careful manual investigation.
    if not previously_archived or \
          (previously_archived and (self.force_overwrite or self.force_md5_check)):

      # Check that the transfer to the archive completed successfully.
      LOGGER.info("Comparing md5 sum of %s in archive and in repository ...",
                  fobj)

      # Raising an exception rolls back the transaction cleanly.
      if not os.path.exists(archpath):
        raise ArchiveError("Error: File has not yet appeared in the archive: %s" % fobj)

      checksum = checksum_file(archpath) # archive path

      if checksum == fobj.checksum:
        LOGGER.info("Md5 sum in repository and for %s are identical.", archpath)

        # Actually record the archiving in the database.
        fobj.save()
        LOGGER.info("Saved file in the %s archive: %s", fobj.archive, fobj)

      else:

        raise ArchiveError(\
          ("Error: Archive file checksum (%s) not same as in"
          + " repository (%s) for file %s. Skipping!") % (checksum, fobj.checksum, fobj))

    return None

  def _remove_older_primary_files_by_type(self):
    '''
    Simple wrapper for _remove_primary_files which identifies primary
    files of the desired filetype and age and queues them for
    deletion.
    '''
    time_threshold = datetime.datetime.now() - \
        datetime.timedelta(days=self.archive.host_delete_timelag)

    if self.filetype == 'fq':
      fobjs = Lanefile.objects.filter(filetype__code=self.filetype,
                                      archive_date__lt=time_threshold)
    elif self.filetype == 'bam':
      fobjs = Alnfile.objects.filter(filetype__code=self.filetype,
                                     archive_date__lt=time_threshold)
    else:
      raise StandardError("Deletion by filetype %s is not supported"
                          % self.filetype) 
    LOGGER.info("Identified %d %s files archived more than %d days ago.",
                len(fobjs), self.filetype, self.archive.host_delete_timelag)

    self._remove_primary_files(fobjs)

  def _remove_primary_files(self, fobjs):
    '''
    Deletes primary copies of the archived files.
    '''
    # If file has been in Archive for long enough or force_delete,
    # delete the source.
    files_deleted = 0
    for fobj in fobjs:

      archpath = fobj.repository_file_path
      repopath = fobj.original_repository_file_path

      if self.force_delete:
        LOGGER.warning(\
          "Executing forced deletion. Archive information: date=%s file=%s."
          + " Removing %s", fobj.archive_date, archpath, repopath)
      else:
        LOGGER.info(\
          "More than %d days passed since archiving %s."
          + " (Archive date=%s).",
          self.archive.host_delete_timelag, repopath, fobj.archive_date)

      # Before deleting the file, check that the file in archive not
      # only exists but has the same md5 sum as recorded in
      # repository. BEAR IN MIND that the way the queries are
      # currently set up this code is run over ALMOST EVERY BAM FILE
      # IN THE REPOSITORY!
      if os.path.exists(repopath):
        if os.path.exists(archpath):
          checksum = checksum_file(archpath)
          if checksum != fobj.checksum:
            LOGGER.error(\
              "Error: Archive file checksum (%s) not same as in repository (%s)."
              + " Can not delete the file!", checksum, fobj.checksum)
          else:
            LOGGER.warning("Removing %s.", repopath)
            os.unlink(repopath)
            files_deleted += 1
        else:
          LOGGER.error("File %s recorded to be in archive but missing on disk.", archpath)
          continue

    if files_deleted == 0:
      LOGGER.info("Zero files deleted from repository.")
    else:
      LOGGER.warning("%d files removed from repository.", files_deleted)

  @transaction.atomic
  def restore_file_from_archive(self, fpath):
    '''
    Method restores the file in the archive back into the main
    repository file tree, and resets the archive flag such that the
    repository copy is now considered authoritative. Note that using
    this function may be counterproductive if running a regular cron
    job to archive all files of a specific filetype (subsequent cron
    job runs will simply move the file back into the archive
    again). It is recommended to only restore files which are not
    managed in this fashion.
    '''
    fobj = _find_file(fpath)

    if not fobj.archive:
      LOGGER.error(\
        "File %s is not currently registered to any archive location.", fobj)

    archpath = fobj.repository_file_path

    # This removes all archive metadata (the repository copy will be
    # authoritative once more). Note that we are not deleting the
    # archived file (since the archive will typically not allow this).
    fobj.archive      = None
    fobj.archive_date = None
    fobj.save()
    repopath = fobj.repository_file_path

    # This may happen when restoring a file halfway through the
    # transfer process. We default to using the copy already present,
    # not least because if e.g. the archive is unavailable this is the
    # only way to restore access quickly.
    if os.path.exists(repopath):
      LOGGER.warning("File %s already present in repository tree.", fobj)
    else:
      copy2(archpath, repopath)

    checksum = checksum_file(repopath)

    # Another manual investigation type error.
    if checksum != fobj.checksum:
      raise ValueError(\
        "Restored file checksum (%s) does not agree with repository value (%s)."
        % (checksum, fobj.checksum))

  def run_archival(self, fnames):
    '''
    Manager method to run bulk archival of files; typically called
    from a cronjob managing the staged transfer of files to the
    archive followed by eventual deletion from the original repository
    (after a suitable time period to allow for archive backup).
    '''
    # In case file type was provided, override the list of files (if
    # any) provided as arguments
    if self.filetype:
      if len(fnames) > 0:
        LOGGER.warning(\
          "Ignoring provided filenames in favour of files of type %s",
          self.filetype)
      fobjs = _get_files_for_filetype(self.filetype, not_archived=True)

      # If we've imposed a lag before archiving (to prevent erroneous
      # storage of the wrong files on a read-only filesystem) we
      # filter files by repository insertion date here.
      if self.archive_lag is not None:
        time_threshold = datetime.datetime.now() - \
            datetime.timedelta(days=self.archive_lag)
        fobjs = fobjs.filter(date__lte=time_threshold)

      LOGGER.info("Found %d non-archived files for copying.", fobjs.count())
    else:

      # This is a little unweildy but since the list may contain a mix
      # of Lanefile, Alnfile and whatever, we cannot simply query
      # Datafile directly.
      fobjs = [ _find_file(fname) for fname in fnames ]

    # From here on, fobjs could either be an array or QuerySet.
    if len(fobjs) > 0:

      # Copy files to the archive. If we plan to wait for the archive
      # filesystem to catch up with reality (due to latency in the
      # system), we first copy everything across as a batch job before
      # we even think about touching the database.
      for fobj in fobjs:
        LOGGER.info("Copying \'%s\' to archive.", fobj)
        self._copy_file_to_archive_disk(fobj)

      # Wait for files copied to archive to become visible in the
      # file system.
      if self.copy_wait_archive:
        LOGGER.info("Copying finished. Waiting 5 minutes for"
                      + " file system to register copied files.")
        time.sleep(5*60)

      # Check files copied successfully, and enter archive information
      # in the database.
      if not self.copy_only:
        LOGGER.warning("Archiving %d non-archived files:", len(fobjs))
        failedfns = []
        for fobj in fobjs:
          try:
            self._register_file_in_archive(fobj)
          except ArchiveError, err:
            LOGGER.error("Archiving failed: %s", err)
            failedfns.append(fobj)

        if len(failedfns) > 0:
          LOGGER.error("%d failed archiving (see above for reasons):",
                       len(failedfns))
          for failedfn in failedfns:
            LOGGER.error("  %s", failedfn)

    # Check for deletion and delete primary copies of files archived
    # long time ago.
    if self.filetype:
      self._remove_older_primary_files_by_type()
    elif self.force_delete:
      LOGGER.warning(\
        "Attention: FORCED DELETION of primary copies for %d archived files!",
        len(fobjs))
      self._remove_primary_files(fobjs)

################################################################################
