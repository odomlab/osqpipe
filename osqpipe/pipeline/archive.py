#!/usr/bin/env python

'''
Code used to manage the transfer and restoration of files to and from
various ArchiveLocations.
'''

import os
import sys
import time
import datetime

from subprocess import Popen, PIPE
from shutil import copy2

from django.db import transaction
from osqpipe.models import ArchiveLocation, Lanefile, Alnfile, \
    QCfile, Peakfile, MergedAlnfile, Datafile
from osqpipe.pipeline.utilities import checksum_file, bash_quote

from osqpipe.pipeline.config import Config
from osqpipe.pipeline.setup_logs import configure_logging

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

  # We will need to double-quote the destination passed to scp (FIXME
  # test this).
  host_archdir = os.path.join(arch.host_path, fobj.libcode)
  dest = bash_quote(os.path.join(host_archdir,
                                 os.path.basename(fobj.repository_file_path)))
  cmd = [ 'scp', '-p',
          '-o', 'StrictHostKeyChecking=no' ]
  if arch.host_port is not None:
    cmd += [ '-P', str(arch.host_port) ]

  # Assume we're copying from the main repository to the archive.
  cmd += [ fobj.original_repository_file_path ]

  if arch.host_user is not None:
    cmd += [ '%s@%s:%s' % (arch.host_user, arch.host, bash_quote(dest)) ]
  else:
    cmd += [ '%s:%s' % (arch.host, bash_quote(dest)) ]

  start_time = time.time()

  while attempts > 0:
    subproc = Popen(cmd, stdout=PIPE, stderr=PIPE)
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
  cmd += [ 'mkdir', bash_quote(folder) ]

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

def _copy_file_to_archive(fobj, fname):
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

    LOGGER.info("Copying %s to the archive.", fname)

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
          fobj = Peakfile.objects.get(filename=fname)
        except Peakfile.DoesNotExist:
          try:
            fobj = MergedAlnfile.objects.get(filename=fname)
          except MergedAlnfile.DoesNotExist:
            raise StandardError(\
              "Datafile %s not found in repository." % fname)

  return fobj

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
  '''

  __slots__ = ('filetype', 'archive', 'copy_only', 'copy_wait_archive',
               'force_delete', 'force_md5_check', 'force_overwrite')

  def __init__(self, filetype=None, archive=CONFIG.default_archive,
               copy_only=False, copy_wait_archive=True,
               force_delete=False, force_md5_check=False,
               force_overwrite=False):
    self.filetype          = filetype
    self.archive           = ArchiveLocation.objects.get(name=archive)
    self.copy_only         = copy_only
    self.copy_wait_archive = copy_wait_archive
    self.force_delete      = force_delete
    self.force_md5_check   = force_md5_check
    self.force_overwrite   = force_overwrite

  ## FIXME this whole method needs splitting into smaller parts so we
  ## don't have (a) a ridiculously long method, and (b) and unwieldy
  ## transaction.
  @transaction.commit_on_success
  def move_file_to_archive(self, fpath):
    '''
    Given a file name (or file path), and the name of an archive as recorded
    in the repository, make sure there is a valid copy of the file in the
    archive and delete the file from the primary repository file tree.
    '''
    fname = os.path.basename(fpath)
    parts = os.path.splitext(fname)
    if parts[1] == CONFIG.gzsuffix:
      fname = parts[0]
    fobj = _find_file(fname)

    # Check if file in archive (both on db and on disk).
    already_in_archive = False
    if fobj.archive:
      already_in_archive = True
      LOGGER.info("File %s already in archive. Date of archiving: %s.",
                  fname, fobj.archive_date)
      if self.force_overwrite:
        fobj.archive = self.archive
        fobj.archive_date = time.strftime('%Y-%m-%d')
        LOGGER.warning("Force overwrite. Updating archive record for %s.",
                       fname)
    else:
      fobj.archive = self.archive
      if not self.copy_only:
        fobj.archive_date = time.strftime('%Y-%m-%d')
        LOGGER.info("Creating archive record for %s.", fname)
      else:
        LOGGER.info("Copying %s to archive but not recording in repository.",
                    fname)

    archpath = fobj.repository_file_path

    # Copy file to archive. In case remote host is provided, scp via
    # remote host. Otherwise copy.
    if self.force_overwrite or not os.path.exists(archpath):
      _copy_file_to_archive(fobj, fname)

    elif not already_in_archive:
      LOGGER.info("File %s already in archive. No need to copy.", fname)

    if self.copy_only:
      return None

    # Errors here will typically need careful manual investigation.
    if (already_in_archive and (self.force_overwrite or self.force_md5_check)) \
          or not already_in_archive:
      LOGGER.info("Comparing md5 sum of %s in archive and in repository ...",
                  fname)
      checksum = checksum_file(archpath)
      if checksum != fobj.checksum:
        LOGGER.error(\
          "Error: Archive file checksum (%s) not same as in"
          + " repository (%s). Skipping!", checksum, fobj.checksum)
        return fname
      else:
        fobj.save() # Do we actually need to save it here, or can we
                    # make the transaction scope smaller FIXME?
      LOGGER.info("Md5 sum in repository and for %s are identical.", archpath)
    return None

  # FIXME this method needs reviewing also.
  def remove_primary_files(self, files=None):
    '''
    Deletes primary copies of the archived files in case forced or
    archived more than specified number of days ago
    (self.archive.host_delete_timelag).  In case list of files is empty,
    identifies all files that have been in archive more than specified
    number of days for a particular file type.
    '''
    # In case empty list of files we look to filetype for guidance.
    if files is None or len(files) == 0:
      time_threshold = datetime.datetime.now() - \
          datetime.timedelta(days=self.archive.host_delete_timelag)
      t_date = datetime.date.today()
      if self.filetype == 'fq':
        files = Lanefile.objects.filter(filetype__code=self.filetype,
                                        archive_date__lt=time_threshold)
      if self.filetype == 'bam':
        files = Alnfile.objects.filter(filetype__code=self.filetype,
                                       archive_date__lt=time_threshold)
      LOGGER.info("Identified %d %s files archived more than %d days ago.",
                  len(files), self.filetype, self.archive.host_delete_timelag)
    if self.force_delete:
      LOGGER.warning(\
        "Attention: FORCED DELETION of primary copies for %d archived files!",
        len(files))

    # If file has been in Archive for long enough or force_delete,
    # delete the source.
    files_deleted = 0
    for fobj in files:

      # Requerying when we already have a Datafile object is just asinine.
      if not issubclass(type(fobj), Datafile):
        fobj = _find_file(fobj) # Assume a str/unicode filename

      archpath = fobj.repository_file_path
      repopath = fobj.original_repository_file_path

      # Convoluted logging code here, should be placed elsewhere
      # rather than confusing the flow with another conditional block. FIXME.
      if self.force_delete:
        LOGGER.warning(\
          "Executing forced deletion. Archive information: date=%s file=%s."
          + " Removing %s", fobj.archive_date, archpath, repopath)
      else:
        LOGGER.info(\
          "More than %d days passed since archiving %s."
          + " (Archive date=%s\tToday=%s).",
          self.archive.host_delete_timelag, repopath, fobj.archive_date, t_date)
        if os.path.exists(repopath):
          LOGGER.warning("Removing %s.", repopath)

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
            os.unlink(repopath)
            files_deleted += 1
      else:
        raise ValueError(\
          "Error: File %s recorded to be in archive but missing on disk."
          % archpath)
    if files_deleted == 0:
      LOGGER.info("No files to delete from repository.")
    else:
      LOGGER.warning("%d files removed from repository.", files_deleted)

  @transaction.commit_on_success
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
    fname = os.path.basename(fpath)
    parts = os.path.splitext(fname)
    if parts[1] == CONFIG.gzsuffix:
      fname = parts[0]
    fobj = _find_file(fname)

    if not fobj.archive:
      LOGGER.error(\
        "File %s is not currently registered to any archive location.", fname)

    archpath = fobj.repository_file_path

    # This removes all archive metadata (the repository copy will be
    # authoritative once more). Note that we are not deleting the
    # archived file (since the archive will typically not allow this).
    fobj.archive      = None
    fobj.archive_date = None
    fobj.save()
    repopath = fobj.repository_file_path

    # This really shouldn't happen, and indicates manual intervention
    # may be necessary.
    if os.path.exists(repopath):
      raise StandardError("File %s already present in repository tree." % fname)

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
      fnames = _get_files_for_filetype(self.filetype, not_archived=True)
      LOGGER.info("Found %d non-archived files for copying.", len(fnames))

    if len(fnames) > 0:
      # Copy files to archive
      if self.copy_wait_archive or self.copy_only:
        for fname in fnames:
          LOGGER.info("Copying \'%s\' to archive.", fname)
          self.move_file_to_archive(str(fname))

      # Wait for files copied to archive to become visible in the file system
      if self.copy_wait_archive:
        LOGGER.info("Copying finished. Waiting 5 minutes for"
                    + " file system to register copied files.")
        time.sleep(5*60)

      # Check files to archive
      failedfns = []
      if not self.copy_only:
        LOGGER.warning("Archiving %d non-archived files:", len(fnames))
        for fname in fnames:
          LOGGER.warning("""Archiving '%s'.""", fname)
          failedfn = self.move_file_to_archive(str(fname))
          if failedfn is not None:
            failedfns.append(failedfn)
      if failedfns:
        LOGGER.error("%d failed archiving due to md5 check sum differences:",
                     len(failedfns))
        for failedfn in failedfns:
          LOGGER.error(failedfn)

    # Check for deletion and delete primary copies of files archived
    # long time ago.
    if self.filetype:
      fnames = []
    if self.force_delete or len(fnames) == 0:
      self.remove_primary_files(fnames)

################################################################################
