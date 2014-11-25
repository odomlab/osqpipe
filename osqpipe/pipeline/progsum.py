#!/usr/bin/env python

'''
Module relating to managing the versions of the external programs
which we use to process our data.
'''

import sys
import os
import subprocess
import re
import logging
from distutils import spawn
from setup_logs import configure_logging
LOGGER = configure_logging()

class ProgramSummary(object):
  """
  This is a class for collecting and holding summary information of a
  linux command line program.  The summary information includes:
  program name, path and version of the program. The information can
  be entered during init but an attempt will be made to validate its
  correctness.
  """
  def __init__(self, program, path = None, version = None,
               versioncmd = None, debug = None,
               ssh_host = None, ssh_user = None, ssh_path = None, ssh_port=None ):

    self.program = None
    self.path    = None
    self.version = None

    self.debug = debug
    if self.debug:
      LOGGER.setLevel(logging.DEBUG)
    else:
      LOGGER.setLevel(logging.INFO)

    if ssh_host is None:
      predictedversion = self.initialize_local(program, path, versioncmd)
    else:
      predictedversion = self.initialize_remote(program, ssh_host, ssh_user,
                                                ssh_path, ssh_port, versioncmd)

    if version is not None and version != "":
      self.version = version

      # Warn if specified and predicted versions are different.
      if version != predictedversion:
        LOGGER.warn(("Specified version \"%s\" different from"
                     + " program version \"%s\""),
                    version, predictedversion)
        self.version = predictedversion
    else:
      self.version = predictedversion

  def initialize_local(self, program, path=None, versioncmd=None):
    '''
    Fairly strict method used to detect program version on the local
    server.
    '''
    # Initialise path; check if program name contains full path.
    (ppath, pprogram) = os.path.split(program)
    if pprogram is None:
      err = "Program name missing!"
      LOGGER.error(err)
      sys.exit(err)

    if ppath is not None and os.path.isdir(ppath) and os.path.isfile(program):
      self.path = ppath

    elif path is not None and os.path.isdir(path) \
          and os.path.isfile(os.path.join(path, pprogram)):
      self.path = path

    else:
      self.path = self.which(pprogram)

    # Check file.
    testpath = os.path.join(self.path, pprogram)
    if os.path.isfile(testpath):
      self.program = pprogram
    else:
      err = ("Program name \"%s\" not a file, unaccesible or missing!"
             % testpath)
      LOGGER.error(err)
      sys.exit(err)

    # Check version.
    predictedversion = self.get_version(versioncmd)
    return predictedversion

  def initialize_remote(self, program, ssh_host,
                        ssh_user, ssh_path, ssh_port, versioncmd=None):
    '''
    Slightly more relaxed method, used by BwaDesktopJobSubmitter,
    to detect the version of the program available to ssh_user on the
    remote ssh_host.
    '''
    (ppath, pprogram) = os.path.split(program)
    self.path    = ppath
    self.program = pprogram
    predictedversion = self.get_version(versioncmd=versioncmd,
                                        ssh_host=ssh_host,
                                        ssh_user=ssh_user,
                                        ssh_path=ssh_path,
                                        ssh_port=ssh_port)
    return predictedversion

  @staticmethod
  def which(program):
    '''
    Identify the installation directory for the program of interest.
    '''
    # FIXME currently assumes os.environ['PATH'] is to be searched; we
    # need to be able to pass in specific path listings.
    out = spawn.find_executable(program)
    (pdir, _prog) = os.path.split(out)

    if os.path.isdir(pdir):
      return pdir

    else:
      LOGGER.error("Path for program %s can not be found!", program)
      sys.exit("Error: Path for program %s can not be found!" % program)

  def get_version(self, versioncmd = None, ssh_host = None,
                  ssh_user = None, ssh_path = None, ssh_port=None):
    '''
    Core entry method for the class. Figure out the version of the
    specified program which is available, either locally or on an SSH
    server depending on wheth ssh_host is set or not.
    '''
    vstrings = ["--version", "-v", " 2>&1 | grep \'[V|v]ersion\'"]
    vpattern = re.compile(r"v?((?:\d+[._-])*\d+)")

    if ssh_host is None:
      program = os.path.join(self.path, self.program)
    else:
      program = ("ssh -p %s %s@%s \"PATH=%s %s\""
                 % (str(ssh_port), ssh_user, ssh_host, ssh_path, self.program))

    # If a string to tease out version from program is known:
    if versioncmd is not None and versioncmd != "":
      vstrings.insert(-1, versioncmd)

    # Try to tease out program version using each of the methods in
    # list until version can be identified.
    for vstring in vstrings:
      version = None
      proc = subprocess.Popen(program + " " + vstring, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, shell=True)
      (out, err) = proc.communicate()
      vmatch = vpattern.search(out)
      if vmatch:
        version = vmatch.group(1)
      else:
        vmatch = vpattern.search(err)
        if vmatch:
          version = vmatch.group(1)
      if version is not None:
        if re.search(r'\.', version):
          # Pull out trailing -\d+ if the version contains periods. This
          # is in part to allow tally/reaper version strings ("13-274"
          # => "13-274") while not overcomplicating bwa/samtools
          # version strings ("0.1.19-44428cd" => "0.1.19")
          version = re.sub(r'-\d+$', '', version)
        return version
    return None

