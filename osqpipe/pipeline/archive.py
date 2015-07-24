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
from osqpipe.models import ArchiveLocation, Lanefile, Alnfile, QCfile, Peakfile, MergedAlnfile
from osqpipe.pipeline.utilities import checksum_file

from osqpipe.pipeline.config import Config
from osqpipe.pipeline.setup_logs import configure_logging

LOGGER = configure_logging('archive')
CONFIG = Config()

################################################################################
class ArchiveManager(object):

  __slots__ = ()

  def __init__(self):
    pass

  
  def _transfer_over_scp(self, source, dest, port=22, user=None,
                          attempts = 1, sleeptime = 2):
    '''
    A wrapper for scp allowing multiple attempts for the transfer in case
    of recoverable error.
    '''
    unrecoverable = [ 'No such file or directory',
                      'Failed to add the host to the list of known hosts',
                      'Operation not permitted' ]
  
    cmd = [ 'scp', '-p',
            '-o', 'StrictHostKeyChecking=no',
            '-P', str(port), source ]
    if user is None:
      cmd += [ dest ]
    else:
      cmd += [ '%s@%s' % (user, dest) ]
  
    t0 = time.time()
  
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
          """Transfer failed with following error code: "%s"\nTrying again (max %d times)""",
          stderr, attempts)
        time.sleep(sleeptime)
      else:
        break
  
    if retcode !=0:
      raise ValueError("ERROR. Failed to transfer %s. (command=\'%s\')\n" % (source, " ".join(cmd)) )
  
    time_diff = time.time() - t0
    LOGGER.info("Copying to archive (scp) completed in %d seconds.", time_diff)
  
    return retcode
  
  def _create_foreign_dir(self, host, port, user, folder):
    '''
    Create folder in foreign host over ssh.
    '''
    # FIXME remove shell=True and use bash_quote etc.
    cmd = 'ssh -p %s %s@%s \'mkdir %s\'' % (port, user, host, folder)
    
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    (stdout, stderr) = p.communicate()
    retcode = p.wait()
    if retcode != 0:
      if 'File exists' in stderr:      
        LOGGER.info("Directory %s@%s:%s already exists.", user, host, folder)
        retcode = 0
      else:
        raise ValueError(\
          "ERROR. Failed to create directory in archive (cmd=\"%s\").\nSTDOUT: %s\nSTDERR: %s\n"
          % (cmd, stdout, stderr) )
  
    return retcode
    
  
  def get_files_for_filetype(self, filetype, not_archived=False):
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
        "'%s' files not supported. Use one of the following file types: [fq, bam]" % filetype)
  
    return files
   
  def _find_file(self, fname): # FIXME reuse code from restoreFileFromArchive
  
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
              raise StandardError("Datafile %s not found in repository." % fname)
  
    return fobj
  
  @transaction.commit_on_success
  def move_file_to_archive(self, fpath, archive, force_overwrite=False,
                            force_delete=False, force_md5_check=False,
                            copy_only=False):
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
    fobj = self._find_file(fname)
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
      LOGGER.info("File %s already in archive. Date of archiving: %s.", fname, fobj.archive_date)
      if force_overwrite:
        fobj.archive = archloc
        fobj.archive_date = time.strftime('%Y-%m-%d')
        # fobj.save() # Do we actually need to save it here, or can we
        # make the transaction scope smaller FIXME?
        LOGGER.warning("Force overwrite. Updating archive record for %s.", fname)
    else:
      fobj.archive = archloc
      if not copy_only:
        fobj.archive_date = time.strftime('%Y-%m-%d') # NB! date format not tested!
        # fobj.save() # Do we actually need to save it here, or can we
        # make the transaction scope smaller FIXME?
        LOGGER.info("Creating archive record for %s.", fname)
      else:
        LOGGER.info("Copying %s to archive but not recording in repository.", fname)
  
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
        retcode = self._create_foreign_dir(archloc.host, archloc.host_port,
                                      archloc.host_user, host_archdir)
        if retcode != 0:
          raise ValueError("Error: Failed to create %s in archive."
                           % (host_archdir))  
        LOGGER.info("Copying %s to the archive.", fname)
        retcode = self._transfer_over_scp(repopath, '%s:%s' % (archloc.host, host_archdir),
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
        copy2(repopath, archpath)
        time_diff = time.time() - t0
        LOGGER.info("Copying to archive completed in %d seconds.", time_diff)
        
    elif not alreadyInArchive:
      LOGGER.info("File %s already in archive. No need to copy.", fname)
    
    if copy_only:
      return None
  
    # Errors here will typically need careful manual investigation.
    if (alreadyInArchive and (force_overwrite or force_md5_check)) or not alreadyInArchive:
      LOGGER.info("Comparing md5 sum of %s in archive and in repository ...", fname)
      checksum = checksum_file(archpath)
      if checksum != fobj.checksum:      
        # raise ValueError("Error: Archive file checksum (%s) not same as in repository (%s)."
        #                  % (checksum, fobj.checksum))
        LOGGER.error("Error: Archive file checksum (%s) not same as in repository (%s). Skipping!",
                     checksum, fobj.checksum)
        return fname
      else:
        fobj.save() # Do we actually need to save it here, or can we make the transaction scope smaller FIXME?
      LOGGER.info("Md5 sum in repository and for %s are identical.", archpath)
    return None
  
  def remove_primary_files(self, files, archive, filetype, force_delete=False):
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
      LOGGER.info("Identified %d %s files archived more than %d days ago.",
                  len(files), filetype, archloc.host_delete_timelag)
    if force_delete:
      LOGGER.warning("Attention: FORCED DELETION of primary copies for %d archived files!",
                     len(files))
  
    # if file has been in Archive for long enough or force_delete, delete the source
    filesDeleted = 0
    for f in files:
      fobj = self._find_file(f)
      archpath = fobj.repository_file_path
      fobj.archive = None
      repopath = fobj.repository_file_path
      if not force_delete:
        LOGGER.info("More than %d days passed since archiving %s. (Archive date=%s\tToday=%s).",
                    archloc.host_delete_timelag, repopath, fobj.archive_date, t_date)
        if os.path.exists(repopath):
          LOGGER.warning("Removing %s.", repopath)
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
      LOGGER.info("No files to delete from repository.")
    else:
      LOGGER.warning("%d files removed from repository.", filesDeleted)
  
  @transaction.commit_on_success
  def restore_file_from_archive(self, fpath):
  
    fname = os.path.basename(fpath)
    parts = os.path.splitext(fname)
    if parts[1] == CONFIG.gzsuffix:
      fname = parts[0]
    fobj = self._find_file(fname)
  
    if not fobj.archive:
      LOGGER.error("File %s is not currently registered to any archive location.", fname)
  
    archpath = fobj.repository_file_path
    fobj.archive = None
    fobj.save()  # FIXME is this necessary? see moveFileToArchive
    repopath = fobj.repository_file_path
  
    # This really shouldn't happen, and indicates manual intervention may be necessary.
    if os.path.exists(repopath):
      raise StandardError("File %s already present in repository tree." % fname)
  
    copy2(archpath, repopath)
  
    checksum = checksum_file(repopath)
  
    # Another manual investigation type error.
    if checksum != fobj.checksum:
      raise ValueError("Restored file checksum (%s) does not agree with repository value (%s)."
                        % (checksum, fobj.checksum))
  
  
################################################################################
