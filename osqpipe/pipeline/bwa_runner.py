#!/usr/bin/env python
#
# $Id$

'''Parent class for the code running BWA on the cluster. Here we
handle some core external dependencies, logging and group
permissions.'''

import sys
import os
import re
import stat
import logging
import grp
import glob
from pipes import quote
from shutil import move
from distutils import spawn
from socket import getfqdn, socket, AF_UNIX, SOCK_STREAM
from tempfile import gettempdir, NamedTemporaryFile
from getpass import getuser

from utilities import call_subprocess, bash_quote, \
    is_zipped, set_file_permissions, BamPostProcessor, parse_repository_filename
from config import Config

from setup_logs import configure_logging
LOGGER = configure_logging('bwa_runner')

##############################################################################

def make_bam_name_without_extension(fqname):

  """Creates bam file basename out of Odom/Carroll lab standard fq
  filename. Note that this function does not append '.bam' to the
  returned name, so that samtools can add it for us."""

  base         = os.path.splitext(fqname)[0]
  lane_pattern = re.compile(r'^(.*)p[12](@\d+)?$')
  matchobj     = lane_pattern.match(base)
  if matchobj != None:
    base = matchobj.group(1)
    if matchobj.group(2):
      base += matchobj.group(2)
  return base

def paired_sanity_check(filenames, is_paired):

  '''Quick argument checking for sanity. We could just count the
  files and set is_paired ourselves, but this allows for a double
  check.'''

  if len(filenames) == 1:
    if is_paired:
      raise ValueError(
        "Only one filename passed for paired-end sequencing alignment.")
  elif len(filenames) == 2:
    if not is_paired:
      raise ValueError(
        "Two filenames passed for single-end sequencing alignment.")
  else:
    raise ValueError(
      "Alignment must be passed either one or two filenames.")

##############################################################################
##############################################################################

class SimpleCommand(object):
  '''
  Simple class used as a default command-string builder.
  '''
  __slots__ = ('conf')

  def __init__(self):

    self.conf = Config()

  def build(self, cmd, *args, **kwargs):
    
    if type(cmd) in (str, unicode):
      cmd = [cmd]

    return " ".join(cmd)

class NohupCommand(SimpleCommand):
  '''
  Class used to wrap a command in a nohup nice invocation for
  execution on a remote desktop.
  '''
  def build(self, cmd, remote_wdir, *args, **kwargs):

    cmd = super(NohupCommand, self).build(cmd, *args, **kwargs)

    # Now we set off a process on the target server which will align
    # the file(s) and scp the result back to us. Note the use of
    # 'nohup' here. Also 'nice', since we're probably using someone's
    # desktop computer. Command is split into two parts for clarity:
    cmd = ("nohup nice -n 20 sh -c '( (%s" % cmd)

    # Closing out the nohup subshell. See this thread for discussion:
    # http://stackoverflow.com/a/29172
    cmd += (") &)' >> %s 2>&1 < /dev/null"
            % os.path.join(remote_wdir, 'remote_worker.log'))

    return cmd

class BsubCommand(SimpleCommand):
  '''
  Class used to build a bsub-wrapped command.
  '''
  def build(self, cmd, mem=2000, queue=None, jobname=None,
            auto_requeue=False, depend_jobs=None, sleep=0, 
            mincpus=1, maxcpus=1, clusterlogdir=None, environ=None, *args, **kwargs):

    # The environ argument allows the caller to pass in arbitrary
    # environmental variables (e.g., JAVA_HOME) as a dict.
    if environ is None:
      environ = {}

    # Pass the PYTHONPATH to the cluster process. This allows us to
    # isolate e.g. a testing instance of the code from production.
    # Note that we can't do this as easily for PATH itself because
    # bsub itself is in a custom location on the cluster.
    for varname in ('PYTHONPATH', 'OSQPIPE_CONFDIR'):
      if varname in os.environ:
        environ[varname] = os.environ[varname]

    cmd = super(BsubCommand, self).build(cmd, *args, **kwargs)

    # Note that if this gets stuck in an infinite loop you will need
    # to use "bkill -r" to kill the job on LSF. N.B. exit code 139 is
    # a core dump. But so are several other exit codes; add 128 to all
    # the unix signals which result in a dump ("man 7 signal") for a
    # full listing.
    qval = "-Q 'all ~0'" if auto_requeue else ''
    try:
      group = self.conf.clustergroup
      if group != '':
        group = '-G ' + group
    except AttributeError:
      group = ''

    resources = 'rusage[mem=%d]' % mem
    memreq    = ''

    # Sanger cluster (farm3) has a stricter set of requirements which
    # are not supported by our local cluster. Quelle surprise.
    try:
      provider = self.conf.clusterprovider
      if provider[:3].lower() == 'san':
        resources = ('select[mem>%d] ' % mem) + resources
        memreq    = '-M %d' % mem
    except AttributeError:
      pass

    # A safety net in case min or max nr of cores gets muddled up. An
    # explicit error is preferred in such cases, so that we can see
    # what to fix.
    if mincpus > maxcpus:
      raise ValueError("mincpus (%d) is greater than maxcpus (%d). Surely some error?" % (mincpus, maxcpus))

    # In case clusterlogdir has been specified, override the self.conf.clusterstdout
    # This is handy in case we want to keep the logs together with job / larger project related files.
    cluster_stdout_stderr = ""
    if clusterlogdir is not None:
      cluster_stdout_stderr = "-o %s/%%J.stdout -e %s/%%J.stderr" % (clusterlogdir, clusterlogdir)
    else:
      cluster_stdout_stderr = "-o %s/%%J.stdout -e %s/%%J.stderr" % (self.conf.clusterstdoutdir, self.conf.clusterstdoutdir)

    envstr  = " ".join([ "%s=%s" % (key, val) for key, val in environ.iteritems() ])
    bsubcmd = (("%s bsub -R '%s' -R 'span[hosts=1]'"
           + " %s -r %s -n %d,%d"
                + " %s %s")
           % (envstr,
              resources,
              memreq,
              cluster_stdout_stderr,
              mincpus,
              maxcpus,
              qval,
              group))

    if queue is not None:
      bsubcmd += ' -q %s' % queue

    # The jobname attribute is also used to control LSF job array creation.
    if jobname is not None:
      bsubcmd += ' -J %s' % jobname

    if depend_jobs is not None:
      depend = "&&".join([ "ended(%d)" % (x,) for x in depend_jobs ])
      bsubcmd += " -w '%s'" % depend

    if sleep > 0:
      cmd = ('sleep %d && ' % sleep) + cmd

    # To group things in a pipe (allowing e.g. use of '&&'), we use a
    # subshell. Note that we quote the sh -c string once, and
    # everything within that string twice. Commands can be of the following form:
    #
    # "ssh blah@blah.org 'echo stuff | target.txt'"
    # r"ssh blah@blah.org \"echo stuff | target.txt\""
    #
    # I.e., one needs to be careful of python's rather idiosyncratic
    # string quoting rules, and use the r"" form where necessary.
    bsubcmd += r' sh -c "(%s)"' % re.sub(r'"', r'\"', cmd)

    return bsubcmd    

##############################################################################
##############################################################################

class JobRunner(object):
  '''
  Instantiable base class used as a core definition of how the various
  job submitter classes are organised. Each JobRunner subclass has a
  command_builder attribute which is used to create the final command
  string that is executed. The default behaviour of this class is to
  simply run the command using call_subprocess on the current
  host. Various subclasses have been created for extending this to
  execute the command on a remote host, or on an LSF cluster. Note
  that to submit to a local LSF head node, one might use this::

  jr = JobRunner(command_builder=BsubCommand())
  jr.submit_command(cmd, mem=8000, queue='dolab')

  See the ClusterJobSubmitter class for how this has been extended to
  submitting to a remote LSF head node.
  '''
  __slots__ = ('test_mode', 'config', 'command_builder')

  def __init__(self, test_mode=False, command_builder=None, *args, **kwargs):
    self.test_mode = test_mode
    if test_mode:
      LOGGER.setLevel(logging.DEBUG)
    else:
      LOGGER.setLevel(logging.INFO)

    self.config = Config()

    self.command_builder = SimpleCommand() \
        if command_builder is None else command_builder

  def run_command(self, cmd, tmpdir=None, path=None, command_builder=None, *args, **kwargs):

    if command_builder:
      cmd = command_builder.build(cmd, *args, **kwargs)
    else:
      cmd = self.command_builder.build(cmd, *args, **kwargs)

    if path is None:
      path = self.config.hostpath

    if tmpdir is None:
      tmpdir = gettempdir()

    LOGGER.debug(cmd)
    if not self.test_mode:
      return call_subprocess(cmd, shell=True, path=path, tmpdir=tmpdir)
    return None

  def submit_command(self, *args, **kwargs):
    '''
    Submit a remote command to whatever queuing or backgrounding
    mechanism the host server supports. Typically this method should
    be used for the big jobs, and run_command for trivial,
    order-sensitive things like checking for the existence of a genome
    on the server, or uncompressing remote files.
    '''
    return self.run_command(*args, **kwargs)

class JobSubmitter(JobRunner):

  '''Class to run jobs via LSF/bsub on the local host (i.e., when running on the cluster).'''

  def __init__(self, remote_wdir=None, *args, **kwargs):
    self.conf = Config()
    super(JobSubmitter, self).__init__(command_builder=BsubCommand(),
                                       *args, **kwargs)

  def submit_command(self, cmd, *args, **kwargs):
    '''
    Submit a job to run on the cluster. Uses bsub to enter jobs into
    the LSF queuing system. Extra arguments are passed to
    BsubCommand.build(). The return value is the integer LSF job ID.
    '''
    pout = super(JobSubmitter, self).\
        submit_command(cmd,
                       *args, **kwargs)

    # FIXME this could be farmed out to utilities?
    jobid_pattern = re.compile(r"Job\s+<(\d+)>\s+is\s+submitted\s+to")
    for line in pout:
      matchobj = jobid_pattern.search(line)
      if matchobj:
        return int(matchobj.group(1))

    raise ValueError("Unable to parse bsub output for job ID.")

class RemoteJobRunner(JobRunner):
  '''
  Abstract base class holding some common methods used by classes
  controlling alignment job submission to the cluster and to other
  computing resources.
  '''
  remote_host = None
  remote_port = 22  # ssh default
  remote_user = None
  remote_wdir = None
  transfer_host = None
  transfer_wdir = None

  def __init__(self, *args, **kwargs):

    # A little programming-by-contract, as it were.
    if not all( x in self.__dict__.keys()
                for x in ('remote_host', 'remote_user', 'remote_wdir',
                          'transfer_host', 'transfer_wdir')):
      raise StandardError("Remote host information not provided.")
    super(RemoteJobRunner, self).__init__(*args, **kwargs)

  def run_command(self, cmd, wdir=None, path=None, command_builder=None, *args, **kwargs):
    '''
    Method used to run a command *directly* on the remote host. No
    attempt will be made to use any kind of queuing or backgrounding
    mechanism.

    The command call is wrapped in an ssh connection. This command will
    also automatically change to the configured remote working
    directory before executing the command.
    '''
    if command_builder:
      cmd = command_builder.build(cmd, *args, **kwargs)
    else:
      cmd = self.command_builder.build(cmd, *args, **kwargs)

    if wdir is None:
      wdir = self.remote_wdir

    if path is None:
      pathdef = ''
    else:
      if type(path) is list:
        path = ":".join(path)
      pathdef = "PATH=%s" % path

    cmd = ("ssh -p %s %s@%s \"source /etc/profile; cd %s && %s %s\""
           % (str(self.remote_port),
              self.remote_user,
              self.remote_host,
              wdir,
              pathdef,
              re.sub(r'"', r'\"', cmd)))
    LOGGER.debug(cmd)
    if not self.test_mode:
      return call_subprocess(cmd, shell=True, path=self.config.hostpath)
    return None

  def remote_copy_files(self, filenames, destnames=None):
    '''
    Copy a set of files across to the remote working directory.
    '''
    if destnames is None:
      destnames = filenames
    if len(filenames) != len(destnames):
      raise ValueError("If used, the length of the destnames list"
                       + " must equal that of the filenames list.")
    for i in range(0, len(filenames)):
      fromfn = filenames[i]
      destfn = destnames[i]

      destfile = os.path.join(self.transfer_wdir, destfn)
      destfile = bash_quote(destfile)

      # Currently we assume that the same login credentials work for
      # both the cluster and the data transfer host. Note that this
      # needs an appropriate ssh key to be authorised on both the
      # transfer host and the cluster host.
      cmd = " ".join(('scp', '-P', str(self.remote_port),
                      '-p', '-q', bash_quote(fromfn),
                      "%s@%s:%s" % (self.remote_user,
                                    self.transfer_host,
                                    quote(destfile))))

      LOGGER.debug(cmd)
      if not self.test_mode:
        call_subprocess(cmd, shell=True, path=self.config.hostpath)

  def remote_uncompress_file(self, fname):
    '''
    Given a remote filename, run gunzip via ssh pipe to uncompress
    it. Return the filename with .gz suffix removed.
    '''
    # Note that we're assuming that the name extensions reflect the
    # compression status.
    destfile = os.path.join(self.remote_wdir, fname)

##  Note that double-quoting here gives undesired results if the
##  filename contains square brackets. If problems recur with other
##  filenames, consider modifying bash_quote to omit the
##  square-bracket quoting.
#    destfile = bash_quote(destfile)

    # Assumes that gzip is in the executable path on the remote server.
    LOGGER.info("Uncompressing remote file %s", fname)
    cmd = " ".join(('gzip -f -d', quote(destfile)))
    self.run_command(cmd, command_builder=SimpleCommand())

    # Remove the .gz extension.
    return os.path.splitext(fname)[0]

  def transfer_data(self, filenames, destnames=None):
    '''
    Convenience method to copy data files across to the server,
    uncompress if necessary, and return the fixed destination
    filenames.
    '''
    if destnames is None:
      destnames = filenames

    # Copy the files across.
    self.remote_copy_files(filenames, destnames)

    # Next, call remote gunzip on any files which need it.
    uncomp_names = []
    for num in range(len(destnames)):

      # We have to test the file we copied over, since we'll be
      # reading its header.
      if is_zipped(filenames[num]):
        uncomp = self.remote_uncompress_file(destnames[num])
      else:
        uncomp = os.path.join(self.remote_wdir, destnames[num])
      uncomp_names.append(uncomp)

    return uncomp_names

  def submit(self, *args, **kwargs):
    '''
    Stub method identifying this as an abstract base class. This is
    typically the primary method which defines the command to be run,
    and which calls self.transfer_data() and
    self.submit_command().
    '''
    raise NotImplementedError(
      "Attempted to submit remote job via an abstract base class.")

##############################################################################

class ClusterJobSubmitter(RemoteJobRunner):

  '''Class to run jobs via LSF/bsub on the cluster.'''

  def __init__(self, remote_wdir=None, *args, **kwargs):

    self.conf        = Config()
    self.remote_host = self.conf.cluster
    self.remote_port = self.conf.clusterport
    self.remote_user = self.conf.clusteruser
    self.remote_wdir = self.conf.clusterworkdir if remote_wdir is None else remote_wdir
    try:
      self.transfer_host = self.conf.transferhost
    except AttributeError, _err:
      LOGGER.debug("Falling back to cluster host for transfer.")
      self.transfer_host = self.remote_host
    try:
      self.transfer_wdir = self.conf.transferdir
    except AttributeError, _err:
      LOGGER.debug("Falling back to cluster remote directory for transfer.")
      self.transfer_wdir = self.remote_wdir

    # Must call this *after* setting the remote host info.
    super(ClusterJobSubmitter, self).__init__(command_builder=BsubCommand(),
                                              *args, **kwargs)

  def submit_command(self, cmd, *args, **kwargs):
    '''
    Submit a job to run on the cluster. Uses bsub to enter jobs into
    the LSF queuing system. Extra arguments are passed to
    BsubCommand.build(). The return value is the integer LSF job ID.
    '''
    pout = super(ClusterJobSubmitter, self).\
        submit_command(cmd,
                       path=self.conf.clusterpath,
                       *args, **kwargs)

    jobid_pattern = re.compile(r"Job\s+<(\d+)>\s+is\s+submitted\s+to")
    if not self.test_mode:
      for line in pout:
        matchobj = jobid_pattern.search(line)
        if matchobj:
          return int(matchobj.group(1))

      raise ValueError("Unable to parse bsub output for job ID.")
    else:
      return 0 # Test mode only.

class ClusterJobRunner(RemoteJobRunner):

  '''Class to run jobs via simple SSH on the cluster.'''

  def __init__(self, remote_wdir=None, *args, **kwargs):

    self.conf        = Config()
    self.remote_host = self.conf.cluster
    self.remote_port = self.conf.clusterport
    self.remote_user = self.conf.clusteruser
    self.remote_wdir = self.conf.clusterworkdir if remote_wdir is None else remote_wdir
    try:
      self.transfer_host = self.conf.transferhost
    except AttributeError, _err:
      LOGGER.debug("Falling back to cluster host for transfer.")
      self.transfer_host = self.remote_host
    try:
      self.transfer_wdir = self.conf.transferdir
    except AttributeError, _err:
      LOGGER.debug("Falling back to cluster remote directory for transfer.")
      self.transfer_wdir = self.remote_wdir

    # Must call this *after* setting the remote host info.
    super(ClusterJobRunner, self).__init__(command_builder=SimpleCommand(),
                                           *args, **kwargs)

##############################################################################

class DesktopJobSubmitter(RemoteJobRunner):
  '''
  Class to run jobs on an alternative alignment host (typically a
  desktop computer with multiple cores). The job is run under nohup
  and nice -20, and a log file is created in the working directory on
  the remote host. The return value is a STDOUT filehandle.
  '''
  def __init__(self, *args, **kwargs):

    self.conf        = Config()
    self.remote_host = self.conf.althost
    self.remote_port = self.conf.althostport
    self.remote_user = self.conf.althostuser
    self.remote_wdir = self.conf.althostworkdir
    self.transfer_host = self.remote_host
    self.transfer_dir  = self.remote_wdir

    # Must call this *after* setting the remote host info.
    super(DesktopJobSubmitter, self).__init__(command_builder=NohupCommand(),
                                              *args, **kwargs)

  def submit_command(self, cmd, *args, **kwargs):
    '''
    Submit a command to run on the designated host desktop
    machine.
    '''
    return super(DesktopJobSubmitter, self).\
        submit_command(cmd,
                       path=self.conf.althostpath,
                       remote_wdir=self.remote_wdir,
                       *args, **kwargs)

##############################################################################
##############################################################################

class AlignmentJobRunner(object):
  '''
  An abstract class from which all alignment job runners inherit. All
  subclasses must instantiate a JobRunner class in the 'job' slot
  prior to calling the base __init__ method.
  '''
  __slots__ = ('finaldir', 'genome', 'job', 'conf', 'samplename')

  job = None
  
  def __init__(self, genome, finaldir='.', samplename=None, *args, **kwargs):

    # A little programming-by-contract, as it were.
#    if not all( hasattr(self, x) for x in ('job')):
#      raise StandardError("JobRunner instance not set.")

    self.conf = Config()

    # Support relative paths as input.
    self.finaldir = os.path.realpath(finaldir)

    # Check if genome exists.
    LOGGER.info("Checking if specified genome file exists.")
    cmd = ("if [ -f %s ]; then echo yes; else echo no; fi" % genome)
    LOGGER.debug(cmd)

    if not self.job.test_mode:
      runjob = ClusterJobRunner(test_mode=self.job.test_mode)
      cmdstdoutfile = runjob.run_command(cmd)
      first_line = cmdstdoutfile.readline()
      first_line = first_line.rstrip('\n')
      if first_line != 'yes':
        raise ValueError("Genome %s unacessible or missing." % genome)

    self.genome = genome
    self.samplename = samplename

  @classmethod
  def genome_path(cls, genome, indexdir, genomedir):
    '''
    Returns the expected path to the fasta file for a given genome index.
    '''
    sciname = genome.scientific_name
    sciname = sciname.replace(" ", "_")
    sciname = sciname.lower()
    gpath   = os.path.join(genomedir, sciname,
                           genome.code, indexdir, "%s.fa" % genome.code)

    return gpath

  @classmethod
  def build_genome_index_path(cls, *args, **kwargs):
    '''Stub method identifying this as an abstract base class.'''
    raise NotImplementedError(
      "Attempted to build genome index path in abstract base class.")

##############################################################################

class BwaClusterJobSubmitter(AlignmentJobRunner):

  '''Class representing the submission of a bwa job to the
  cluster. This class in fact uploads the fastq file, gunzips it if
  necessary, and then submits a job to split the fastq file into
  chunks and run the alignments as secondary jobs, also spawning one
  last job which waits for the first to complete before merging the
  output and copying it back to the source server.'''

  def __init__(self, *args, **kwargs):
    self.job = ClusterJobSubmitter(*args, **kwargs)
    super(BwaClusterJobSubmitter, self).__init__(*args, **kwargs)

  def submit(self, filenames,
             is_paired=False, destnames=None, cleanup=True,
             nocc=None, *args, **kwargs):

    '''Actually submit the job. The optional destnames argument can be
    used to name files on the cluster differently to the source. This
    is occasionally useful.'''

    paired_sanity_check(filenames, is_paired)

    # First, copy the files across and uncompress on the server.
    LOGGER.info("Copying files to the cluster.")
    destnames = self.job.transfer_data(filenames, destnames)

    # Next, create flag for cleanup
    if cleanup:
      cleanupflag = '--cleanup'
    else:
      cleanupflag = ''

    # Next, create flag for number of non-unique reads to keep in samse/sampe
    if nocc:
      noccflag = '--n_occ %s' % (nocc,)
    else:
      noccflag = ''

    if self.samplename:

      # Sample names containing spaces are bad on the command line,
      # and potentially problematic in bam read groups.
      sanity_re = re.compile(r'([ \/\(\);&|]+)')
      sampleflag = '--sample %s' % (sanity_re.sub('_', self.samplename),)
    else:
      sampleflag = ''

    # FIXME assumes path on localhost is same as path on cluster.
    progpath = spawn.find_executable('cs_runBwaWithSplit.py', path=self.conf.clusterpath)

    # Next, submit the actual jobs on the actual cluster.
    if is_paired:
      LOGGER.debug("Running bwa on paired-end sequencing input.")
      fnlist = " ".join([ quote(x) for x in destnames ])
      ## FIXME think about ways this could be improved.
      ## In the submitted command:
      ##   --rcp       is where cs_runBwaWithSplit_Merge.py eventually copies
      ##                 the reassembled bam file (via scp).
      cmd = ("python %s --loglevel %d %s %s --rcp %s:%s %s %s %s"
             % (progpath,
                LOGGER.getEffectiveLevel(),
                cleanupflag,
                noccflag,
                self.conf.datahost,
                self.finaldir,
                sampleflag,
                self.genome,
                fnlist))

    else:
      LOGGER.debug("Running bwa on single-end sequencing input.")
      fnlist = quote(destnames[0])
      cmd = ("python %s --loglevel %d %s %s --rcp %s:%s %s %s %s"
             % (progpath,
                LOGGER.getEffectiveLevel(),
                cleanupflag,
                noccflag,
                self.conf.datahost,
                self.finaldir,
                sampleflag,
                self.genome,
                fnlist))

    LOGGER.info("Submitting bwa job to cluster.")
    self.job.submit_command(cmd, *args, **kwargs)

  @classmethod
  def build_genome_index_path(cls, genome, *args, **kwargs):

    conf = Config()

    from progsum import ProgramSummary
    from ..models import Program

    # Get information about default aligner, check that the program is
    # in path and try to predict its version.
    alignerinfo = ProgramSummary(conf.aligner,
                                 ssh_host=conf.cluster,
                                 ssh_port=conf.clusterport,
                                 ssh_user=conf.clusteruser,
                                 ssh_path=conf.clusterpath)
    indexdir = None

    # Check that the version of aligner has been registered in
    # repository.
    try:
      Program.objects.get(program=alignerinfo.program,
                          version=alignerinfo.version,
                          current=True)
      indexdir = "%s-%s" % (alignerinfo.program, alignerinfo.version)

    # If aligner version is missing, try to insert it into the database
    # (FIXME not yet implemented while we see how this works).
    except Program.DoesNotExist, _err:
      sys.exit(("""Aligner "%s" version "%s" found at path "%s" """
               % (alignerinfo.program, alignerinfo.version, alignerinfo.path))
               + "not recorded as current in repository! Quitting.")

    gpath = cls.genome_path(genome, indexdir, conf.clustergenomedir)

    return gpath

##############################################################################

class TophatClusterJobSubmitter(AlignmentJobRunner):
  '''
  Class representing the submission of a tophat2 job to the
  cluster. This class works similarly to BwaClusterJobSubmitter.
  '''
  def __init__(self, *args, **kwargs):
    self.job = ClusterJobSubmitter(*args, **kwargs)
    super(TophatClusterJobSubmitter, self).__init__(*args, **kwargs)

  def submit(self, filenames, is_paired=False, destnames=None, cleanup=True,
              *args, **kwargs):
    '''
    Actually submit the job. The optional destnames argument can be
    used to name files on the cluster differently to the source. This
    is occasionally useful.
    '''
    paired_sanity_check(filenames, is_paired)

    # First, copy the files across and uncompress on the server.
    LOGGER.info("Copying files to the cluster.")
    destnames = self.job.transfer_data(filenames, destnames)

    # Next, create flag for cleanup
    if cleanup:
      cleanupflag = '--cleanup'
    else:
      cleanupflag = ''

    if self.samplename:
      sampleflag = '--sample %s' % (self.samplename,)
    else:
      sampleflag = ''

    # FIXME assumes path on localhost is same as path on cluster.
    progpath = spawn.find_executable('cs_runTophatWithSplit.py', path=self.conf.clusterpath)

    # Next, submit the actual jobs on the actual cluster.
    fnlist = " ".join([ quote(x) for x in destnames ])
    cmd = ("python %s --loglevel %d %s --rcp %s:%s %s %s %s"
           % (progpath,
              LOGGER.getEffectiveLevel(),
              cleanupflag,
              self.conf.datahost,
              self.finaldir,
              sampleflag,
              self.genome,
              fnlist))

    LOGGER.info("Submitting tophat job to cluster.")
    self.job.submit_command(cmd, *args, **kwargs)

  @classmethod
  def build_genome_index_path(cls, genome, *args, **kwargs):

    conf = Config()

    from progsum import ProgramSummary
    from ..models import Program

    # Get information about default aligner, check that the program is
    # in path and try to predict its version.
    alignerinfo = ProgramSummary('bowtie2',
                                 ssh_host=conf.cluster,
                                 ssh_port=conf.clusterport,
                                 ssh_user=conf.clusteruser,
                                 ssh_path=conf.clusterpath)
    indexdir = None

    # Check that the version of aligner has been registered in
    # repository.
    try:
      Program.objects.get(program=alignerinfo.program,
                          version=alignerinfo.version,
                          current=True)
      indexdir = "%s-%s" % ('bowtie', alignerinfo.version)

    except Program.DoesNotExist, _err:
      sys.exit(("""Aligner "%s" version "%s" found at path "%s" """
               % (alignerinfo.program, alignerinfo.version, alignerinfo.path))
               + "not recorded as current in repository! Quitting.")

    # Tophat/bowtie need the trailing .fa removed.
    gpath = cls.genome_path(genome, indexdir, conf.clustergenomedir)

    return gpath

##############################################################################

class BwaDesktopJobSubmitter(AlignmentJobRunner):

  '''An alternative means of submitting an alignment job to a remote
  server. This was originally conceived as a way of working round
  periods of extended downtime on our compute cluster. Jobs may be
  submitted to e.g. a desktop machine with multiple cores and bwa run
  in parallel mode.'''

  def __init__(self, num_threads=1, *args, **kwargs):
    self.job = DesktopJobSubmitter(*args, **kwargs)
    super(BwaDesktopJobSubmitter, self).__init__(*args, **kwargs)
    self.num_threads = num_threads

  def submit(self, filenames,
             is_paired=False, destnames=None, cleanup=True,
             nocc=None, *args, **kwargs):

    '''Submit a job as a background chained process on the designated
    remote server.'''

    paired_sanity_check(filenames, is_paired)

    # First, copy the files across and uncompress on the server.
    LOGGER.info("Copying files to the alignment server.")
    destnames = self.job.transfer_data(filenames, destnames)

    outfnbase = make_bam_name_without_extension(destnames[0])
    outfnfull = outfnbase + '.bam'

    tempfiles = destnames + [ outfnfull ]

    # Next, create flag for number of non-unique reads to keep in samse/sampe
    if nocc:
      if is_paired:
        noccflag = '-o %s' % nocc
      else:
        noccflag = '-n %s' % nocc
    else:
      noccflag = ''

    cmd = ''

    if is_paired:
      LOGGER.debug("Running bwa on paired-end sequencing input.")
      num_threads = max(1, int(self.num_threads))
      fnlist = " ".join([ quote(x) for x in destnames ])

      # We have to run bwa aln on both inputs separately, then combine
      # them with bwa sampe.
      fq_0 = quote(destnames[0])
      fq_1 = quote(destnames[1])

      sai_0 = "%s.sai" % fq_0
      sai_1 = "%s.sai" % fq_1
      tempfiles.extend([ sai_0, sai_1 ])

      # FIXME note that the PATH is not set correctly for pipe sinks
      # in this command (e.g. samtools). Either we need to wrap this
      # somehow (see SplitBwaRunner) or we need to drop PATH in after
      # each '|'.
      cmd += ("bwa aln -t %d %s %s > %s"
              % (num_threads,
                 self.genome,
                 fq_0,
                 sai_0))
      cmd += (" && bwa aln -t %d %s %s > %s"
              % (num_threads,
                 self.genome,
                 fq_1,
                 sai_1))
      cmd += ((" && bwa sampe %s %s %s %s %s %s | samtools view -b -S -u -@ %d"
                                          + " - | samtools sort -@ %d - %s")
              % (noccflag,
                 self.genome,
                 sai_0,
                 sai_1,
                 fq_0,
                 fq_1,
                 num_threads,
                 num_threads,
                 outfnbase))
    else:
      LOGGER.debug("Running bwa on single-end sequencing input.")
      num_threads = max(1, int(self.num_threads))
      fnlist = quote(destnames[0])
      cmd += (("bwa aln -t %d %s %s | bwa samse %s %s - %s"
                                + " | samtools view -b -S -u -@ %d -"
                                + " | samtools sort -@ %d - %s")
             % (num_threads,
                self.genome,
                fnlist,
                noccflag,
                self.genome,
                fnlist,
                num_threads,
                num_threads,
                outfnbase))

    # This is invariant PE vs. SE. First, run our standard picard cleanup:
    postproc = BamPostProcessor(input_fn=outfnfull, output_fn=outfnfull,
                                samplename=self.samplename,
                                tmpdir=self.conf.althostworkdir)
    cmd += (" && %s && rm %s"
            % (" ".join(postproc.clean_sam()), outfnfull))
    cmd += (" && %s && rm %s"
            % (" ".join(postproc.add_or_replace_read_groups()), postproc.cleaned_fn))
    cmd += (" && %s && rm %s"
            % (" ".join(postproc.fix_mate_information()), postproc.rgadded_fn))

    # We generate verbose logging here to better monitor file transfers.
    cmd += (" && scp -v -i %s %s %s@%s:%s"
            % (self.conf.althostsshkey,
               outfnfull,
               self.conf.user,
               self.conf.datahost,
               self.finaldir))
    cmd += (" && ssh -i %s %s@%s touch %s/%s.done"
            % (self.conf.althostsshkey,
               self.conf.user,
               self.conf.datahost,
               self.finaldir,
               outfnfull))

    if cleanup:
      cmd += (" && rm %s" % " ".join(tempfiles))

    LOGGER.info("Submitting bwa job to alignment host.")
    self.job.submit_command(cmd)

  @classmethod
  def build_genome_index_path(cls, genome, *args, **kwargs):

    from progsum import ProgramSummary
    from ..models import Program

    conf = Config()

    # Get information about default aligner, check that the program is
    # in path and try to predict its version.
    alignerinfo = ProgramSummary(conf.aligner,
                                 ssh_host=conf.althost,
                                 ssh_user=conf.althostuser,
                                 ssh_path=conf.althostpath,
                                 ssh_port=conf.althostport)

    # Check that the version of aligner has been registered in
    # repository.
    try:
      Program.objects.get(program=alignerinfo.program,
                          version=alignerinfo.version,
                          current=True)
      indexdir = "%s-%s" % (alignerinfo.program, alignerinfo.version)

    # If aligner version is missing, try to insert it into the database
    # (FIXME not yet implemented while we see how this works).
    except Program.DoesNotExist, _err:
      sys.exit(("""Aligner "%s" version "%s" found at path "%s" """
                % (alignerinfo.program, alignerinfo.version, alignerinfo.path))
               + "not recorded as current in repository! Quitting.")

    gpath = cls.genome_path(genome, indexdir, conf.althostgenomedir)

    return gpath

##############################################################################
##############################################################################

class AlignmentManager(object):
  '''
  Parent class handling various functions required by scripts
  running BWA or other aligners across multiple parallel alignments
  (and merging their output) on the cluster.
  '''
  __slots__ = ('conf', 'samtools_prog', 'group', 'cleanup', 'loglevel',
               'split_read_count', 'bsub', 'merge_prog', 'logfile', 'debug')

  def __init__(self, merge_prog=None, cleanup=False, group=None,
               split_read_count=1000000,
               loglevel=logging.WARNING, debug=True):

    self.conf = Config()

    self.samtools_prog = 'samtools'
    
    # The merge_prog argument *must* be set when calling split_and_align.
    self.merge_prog    = merge_prog
    self.logfile       = self.conf.splitbwarunlog
    self.bsub          = JobSubmitter()
    self.split_read_count   = split_read_count

    self.cleanup       = cleanup
    self.group         = group
    self.debug         = debug
    self._configure_logging(self.__class__.__name__, LOGGER)
    LOGGER.setLevel(loglevel)
    LOGGER.debug("merge_prog set to %s", self.merge_prog)

  def _configure_logging(self, name, logger=LOGGER):
    '''
    Configures the logs to be saved in self.logfile.
    '''
    if self.debug:
      logger.setLevel(logging.DEBUG) # log everything on debug level
    else:
      logger.setLevel(logging.INFO)

    # specify the format of the log file
    logfmt = "[%%(asctime)s] %s %%(levelname)s : %%(message)s" % (name,)
    fmt = logging.Formatter(logfmt)

    # Push stderr to logs; Note that any required StreamHandlers will
    # have been added in the child class.
    hdlr = logging.FileHandler(self.logfile)
    hdlr.setFormatter(fmt)
    hdlr.setLevel(min(logger.getEffectiveLevel(), logging.WARN))
    logger.addHandler(hdlr)

  def split_fq(self, fastq_fn):
    '''
    Splits fastq file to self.split_read_count reads per file using
    linux command line split for speed.
    '''
    LOGGER.debug("splitting fq file %s to %s per file ", fastq_fn, self.split_read_count)

    fastq_fn_suffix = fastq_fn + '-'
    cmd = ("split -l %s %s %s" # split -l size file.fq prefix
           % (self.split_read_count*4, quote(fastq_fn), quote(fastq_fn_suffix)))
    call_subprocess(cmd, shell=True,
                   tmpdir=self.conf.clusterworkdir,
                   path=self.conf.clusterpath)

    # glob will try and expand [, ], ? and *; we don't actually want that.
    # Here we quote them as per the glob docs in a character class [].
    bash_re  = re.compile(r'([?\[\]*])')
    fq_files = glob.glob(bash_re.sub(r'[\1]', fastq_fn_suffix) + "??")
    fq_files.sort()
    for fname in fq_files:
      LOGGER.debug("Created fastq file: '%s'", fname)
      if self.group != None:
        set_file_permissions(self.group, fname)
    if self.cleanup:
      os.unlink(fastq_fn)
      LOGGER.info("Unlinking fq file '%s'", fastq_fn)
    return fq_files

  def queue_merge(self, bam_files, depend, bam_fn, rcp_target, samplename=None):
    '''
    Submits samtools job for merging list of bam files to LSF cluster.
    '''
    assert( self.merge_prog is not None )
    input_files = " ".join(bam_files) # singly-bash-quoted
    LOGGER.debug("Entering queue_merge with input_files=%s", input_files)

    # The self.merge_prog command is assumed to be a python script
    # conforming to the set of arguments supported by
    # cs_runBwaWithSplit_Merge.py. We call 'python' here to pick up
    # the python in our path rather than the python in the merge_prog
    # script shebang line.
    cmd = ("python %s --loglevel %d"
           % (self.merge_prog, LOGGER.getEffectiveLevel()))
    if rcp_target:
      cmd += " --rcp %s" % (rcp_target,)
    if self.cleanup:
      cmd += " --cleanup"
    if self.group:
      cmd += " --group %s" % (self.group,)
    if samplename:
      cmd += " --sample %s" % (samplename,)
    cmd += " %s %s" % (bash_quote(bam_fn), input_files)

    LOGGER.info("preparing samtools merge on '%s'", input_files)
    LOGGER.debug(cmd)

    jobname = bam_files[0].split("_")[0] + "bam"
    jobid = self._submit_lsfjob(cmd, jobname, depend, mem=10000)
    LOGGER.debug("got job id '%s'", jobid)

  def _submit_lsfjob(self, command, jobname, depend=None, sleep=0, mem=8000):
    '''
    Executes command in LSF cluster.
    '''
    jobid = self.bsub.submit_command(command, jobname=jobname,
                                     depend_jobs=depend, mem=mem,
                                     path=self.conf.clusterpath,
                                     tmpdir=self.conf.clusterworkdir,
                                     queue=self.conf.clusterqueue,
                                     sleep=sleep)
    return '' if jobid is None else jobid

  def split_and_align(self, *args, **kwargs):
    '''
    Method used to launch the initial file splitting and bwa
    alignments. This class also submits a job dependent on the outputs
    of those alignments, which in turn merges the outputs to generate
    the final bam file.
    '''
    raise NotImplementedError()

  def _merge_files(self, output_fn, input_fns):
    '''
    Merges list of bam files.
    '''
    if len(input_fns) == 1:
      LOGGER.warn("renaming file: %s", input_fns[0])
      move(input_fns[0], output_fn)
    else:
      cmd = ("%s merge %s %s" # assumes sorted input bams.
             % (self.samtools_prog,
                bash_quote(output_fn),
                " ".join([ bash_quote(x) for x in input_fns])))
      LOGGER.debug(cmd)
      pout = call_subprocess(cmd, shell=True,
                            tmpdir=self.conf.clusterworkdir,
                            path=self.conf.clusterpath)
      for line in pout:
        LOGGER.warn("SAMTOOLS: %s", line[:-1])
    if not os.path.isfile(output_fn):
      LOGGER.error("expected output file '%s' cannot be found.", output_fn)
      sys.exit("File access error.")
    if self.group:
      set_file_permissions(self.group, output_fn)
    if self.cleanup:
      for fname in input_fns:

        # Remove input bam file, as long as it's not the only one.
        if len(input_fns) > 1:
          LOGGER.info("Unlinking bam file '%s'", fname)
          os.unlink(fname)

  def copy_result(self, target, fname):
    '''
    Copies file to target location.
    '''
    qname = bash_quote(fname)
    cmd = "scp -p -q %s %s" % (qname, target)
    LOGGER.debug(cmd)
    pout = call_subprocess(cmd, shell=True,
                          tmpdir=self.conf.clusterworkdir,
                          path=self.conf.clusterpath)
    count = 0
    for line in pout:
      LOGGER.warn("SCP: %s", line[:-1])
      count += 1
    if count > 0:
      LOGGER.error("Got errors from scp, quitting.")
      sys.exit("No files transferred.")
    flds = target.split(":")
    if len(flds) == 2: # there's a machine and path
      fn_base = os.path.basename(qname)
      cmd = "ssh %s touch %s/%s.done" % (flds[0], flds[1], bash_quote(fn_base))
      LOGGER.debug(cmd)
      call_subprocess(cmd, shell=True,
                     tmpdir=self.conf.clusterworkdir,
                     path=self.conf.clusterpath)
    if self.cleanup:
      os.unlink(fname)
    return

  def picard_cleanup(self, output_fn, input_fn, samplename=None):
    '''
    Run picard CleanSam, AddOrReplaceReadGroups,
    FixMateInformation. Note that this method relies on the presence
    of a wrapper shell script named 'picard' in the path.
    '''
    postproc = BamPostProcessor(input_fn=input_fn, output_fn=output_fn,
                                samplename=samplename,
                                tmpdir=self.conf.clusterworkdir)

    # Run CleanSam
    call_subprocess(postproc.clean_sam(),
                    tmpdir=self.conf.clusterworkdir, path=self.conf.clusterpath)
    if self.cleanup:
      os.unlink(input_fn)

    # Run AddOrReplaceReadGroups
    call_subprocess(postproc.add_or_replace_read_groups(),
                    tmpdir=self.conf.clusterworkdir, path=self.conf.clusterpath)
    if self.cleanup:
      os.unlink(postproc.cleaned_fn)
    
    # Run FixMateInformation
    call_subprocess(postproc.fix_mate_information(),
                    tmpdir=self.conf.clusterworkdir, path=self.conf.clusterpath)
    if self.cleanup:
      os.unlink(postproc.rgadded_fn)
    
    if self.group:
      set_file_permissions(self.group, output_fn)

  def merge_alignments(self, input_fns, output_fn, rcp_target=None, samplename=None):
    '''
    Method used to merge a set of bam files into a single output bam
    file.
    '''
    merge_fn = "%s_dirty.bam" % os.path.splitext(output_fn)[0]

    LOGGER.info("merging '%s' into '%s'", ", ".join(input_fns), merge_fn)
    self._merge_files(merge_fn, input_fns)
    LOGGER.info("merged '%s' into '%s'", ", ".join(input_fns), merge_fn)

    LOGGER.info("running picard cleanup on '%s'", merge_fn)
    self.picard_cleanup(output_fn, merge_fn, samplename)
    LOGGER.info("ran picard cleanup on '%s' creating '%s'", merge_fn, output_fn)

    if rcp_target:
      self.copy_result(rcp_target, output_fn)
      LOGGER.info("copied '%s' to '%s'", output_fn, rcp_target)

##############################################################################

class BwaAlignmentManager(AlignmentManager):
  '''
  Subclass of AlignmentManager implementing the bwa-specific
  components of our primary alignment pipeline.
  '''
  def __init__(self, nocc=None, *args, **kwargs):
    super(BwaAlignmentManager, self).__init__(*args, **kwargs)

    # These are now identified by passing in self.conf.clusterpath to
    # the remote command.
    self.bwa_prog      = 'bwa'

    if nocc:
      self.nocc = "-n %s" % (nocc,)
    else:
      self.nocc = ''
    
  def run_bwas(self, genome, paired, fq_files, fq_files2):
    '''
    Submits bwa alignment jobs for list of fq files to LSF cluster.
    '''
    job_ids = []
    out_names = []
    current = 0
    # splits the fq_file by underscore and returns first element which
    # in current name
    for fqname in fq_files:
      donumber = fqname.split("_")[0]
      out = bash_quote(fqname + ".bam")
      out_names.append(out)
      jobname_bam = "%s_%s_bam" % (donumber, current)

      if paired:
        jobname1 = "%s_%s_sai1" % (donumber, current)
        jobname2 = "%s_%s_sai2" % (donumber, current)
        sai_file1 = "%s.sai" % fqname
        sai_file2 = "%s.sai" % fq_files2[current]

        # Run bwa aln
        cmd1 = "%s aln %s %s > %s" % (self.bwa_prog, genome,
                                      bash_quote(fqname),
                                      bash_quote(sai_file1))
        cmd2 = "%s aln %s %s > %s" % (self.bwa_prog, genome,
                                      bash_quote(fq_files2[current]),
                                      bash_quote(sai_file2))

        # Run bwa sampe
        cmd3  = ("%s sampe %s %s %s %s %s %s"
                 % (self.bwa_prog, self.nocc, genome, bash_quote(sai_file1),
                    bash_quote(sai_file2), bash_quote(fqname), bash_quote(fq_files2[current])))

        # Convert to bam
        cmd3 += (" | %s view -b -S -u - > %s.unsorted" % (self.samtools_prog, out))

        # Sort the bam
        cmd3 += (" && %s sort %s.unsorted %s" % (self.samtools_prog, out, bash_quote(fqname)))

        # Cleanup
        cmd3 += (" && rm %s %s %s %s %s.unsorted"
                 % (bash_quote(sai_file1), bash_quote(sai_file2),
                    bash_quote(fqname), bash_quote(fq_files2[current]), out))

        LOGGER.info("starting bwa step1 on '%s'", fqname)
        jobid_sai1 = self._submit_lsfjob(cmd1, jobname1, sleep=current)
        LOGGER.debug("got job id '%s'", jobid_sai1)
        LOGGER.info("starting bwa step1 on '%s'", fq_files2[current])
        jobid_sai2 = self._submit_lsfjob(cmd2, jobname2, sleep=current)
        LOGGER.debug("got job id '%s'", jobid_sai2)

        if jobid_sai1 and jobid_sai2:
          LOGGER.info("preparing bwa step2 on '%s'", fqname)
          jobid_bam = self._submit_lsfjob(cmd3, jobname_bam,
                                          (jobid_sai1, jobid_sai2), sleep=current)
          LOGGER.debug("got job id '%s'", jobid_bam)
          job_ids.append(jobid_bam)
        else:
          LOGGER.error("bjob submission for bwa step1 for '%s' or '%s' failed!",
                       fqname, fq_files2[current])
      else:

        # Run bwa aln
        cmd  = ("%s aln %s %s" % (self.bwa_prog, genome, bash_quote(fqname)))

        # Run bwa samse
        cmd += (" | %s samse %s %s - %s" % (self.bwa_prog, self.nocc,
                                            genome, bash_quote(fqname)))
        # Convert to bam
        cmd += (" | %s view -b -S -u - > %s.unsorted" % (self.samtools_prog, out))

        # Sort the output bam
        cmd += (" && %s sort %s.unsorted %s" % (self.samtools_prog,
                                                out, bash_quote(fqname)))
        # Clean up
        cmd += (" && rm %s %s.unsorted" % (bash_quote(fqname), out))
        
        LOGGER.info("starting bwa on '%s'", fqname)
        LOGGER.debug(cmd)
        jobid_bam = self._submit_lsfjob(cmd, jobname_bam, sleep=current)
        LOGGER.debug("got job id '%s'", jobid_bam)
        job_ids.append(jobid_bam)
      current += 1

    return (job_ids, out_names)

  def split_and_align(self, files, genome, samplename, rcp_target=None):
    '''
    Method used to launch the initial file splitting and bwa
    alignments. This class also submits a job dependent on the outputs
    of those alignments, which in turn merges the outputs to generate
    the final bam file.
    '''
    assert( self.merge_prog is not None )
    fq_files = self.split_fq(files[0])
    paired = False
    if len(files) == 2:
      fq_files2 = self.split_fq(files[1])
      paired = True
    elif len(files) == 1:
      fq_files2 = None
    else:
      LOGGER.error("Too many files specified.")
      sys.exit("Unexpected number of files passed to script.")

    (job_ids, bam_files) = self.run_bwas(genome, paired, fq_files, fq_files2)

    bam_fn = "%s.bam" % make_bam_name_without_extension(files[0])
    self.queue_merge(bam_files, job_ids, bam_fn, rcp_target, samplename)

class TophatAlignmentManager(AlignmentManager):
  '''
  Subclass of AlignmentManager implementing the tophat2-specific
  components of our primary alignment pipeline.
  '''
  def __init__(self, *args, **kwargs):
    super(TophatAlignmentManager, self).__init__(*args, **kwargs)

    # These are now identified by passing in self.conf.clusterpath to
    # the remote command.
    self.tophat_prog   = 'tophat2'
    
  def run_tophat(self, genome, paired, fq_files, fq_files2):
    '''
    Submits tophat2 alignment jobs for list of fq files to LSF cluster.
    '''
    job_ids = []
    out_names = []
    current = 0

    # Tophat/bowtie requires the trailing .fa to be removed.
    genome = re.sub(r'\.fa$', '', genome)

    for fqname in fq_files:
      (donumber, facility, lanenum, _pipe) = parse_repository_filename(fqname)

      # Used as a job ID and also as an output directory, so we want
      # it fairly collision-resistant.
      if donumber is None: # unusual filename, non-repository.
        jobname_bam = "%s_tophat" % fqname
      else:
        jobname_bam = "%s_%s%02d_%s_tophat" % (donumber, facility, int(lanenum), current)

      out = bash_quote(fqname + ".bam")
      out_names.append(out)

      # Run tophat2. The no-coverage-search option is required when
      # splitting the fastq file across multiple cluster nodes. The
      # fr-firststrand library type is the Odom lab default. We
      # use the -p option to ask for more threads; FIXME config option?
      cmd  = ("%s --no-coverage-search --library-type fr-firststrand -p 4 -o %s %s %s"
               % (self.tophat_prog, jobname_bam, genome, bash_quote(fqname)))
      if paired:
        cmd += " %s" % (bash_quote(fq_files2[current]),)

      # Merge the mapped and unmapped outputs, clean out unwanted
      # secondary alignments. Tophat2 sorts the output bams by default.
      strippedbam = "%s.partial" % out
      cmd += (" && %s view -b -F 0x100 -o %s %s"
               % (self.samtools_prog, strippedbam,
                   os.path.join(jobname_bam, 'accepted_hits.bam')))
      cmd += (" && %s merge %s %s %s"
               % (self.samtools_prog, out, strippedbam,
                  os.path.join(jobname_bam, 'unmapped.bam')))

      # Clean up
      cmd += (" && rm -r %s %s %s" % (jobname_bam, strippedbam, bash_quote(fqname)))
      if paired:
        cmd += " %s" % (bash_quote(fq_files2[current]),)
        
      LOGGER.info("starting tophat2 on '%s'", fqname)
      LOGGER.debug(cmd)
      jobid_bam = self._submit_lsfjob(cmd, jobname_bam, sleep=current)
      LOGGER.debug("got job id '%s'", jobid_bam)
      job_ids.append(jobid_bam)

      current += 1

    return (job_ids, out_names)
    
  def split_and_align(self, files, genome, samplename, rcp_target=None):
    '''
    Method used to launch the initial file splitting and bwa
    alignments. This class also submits a job dependent on the outputs
    of those alignments, which in turn merges the outputs to generate
    the final bam file.
    '''
    assert( self.merge_prog is not None )
    fq_files = self.split_fq(files[0])
    paired = False
    if len(files) == 2:
      fq_files2 = self.split_fq(files[1])
      paired = True
    elif len(files) == 1:
      fq_files2 = None
    else:
      LOGGER.error("Too many files specified.")
      sys.exit("Unexpected number of files passed to script.")
      
    (job_ids, bam_files) = self.run_tophat(genome, paired, fq_files, fq_files2)

    bam_fn = "%s.bam" % make_bam_name_without_extension(files[0])
    self.queue_merge(bam_files, job_ids, bam_fn, rcp_target, samplename)


##############################################################################
##############################################################################

class ClusterJobManager(object):
  '''
  Moderately abstract base class providing some methods and attributes
  commonly used in higher-level cluster process management classes
  (e.g. GsnapManager, LastzManager).
  '''
  __slots__ = ('namespace', 'submitter', 'runner', 'config',
               'throttle', 'memsize', 'ssh_key', 'local_workdir')

  def __init__(self, namespace=None, throttle=0, memsize=20,
               ssh_key=None, local_workdir='.'):

    self.config = Config()

    if namespace is None:
      self.namespace = str(os.getpid())
    else:
      self.namespace = namespace

    # These will default to the config cluster working directory.
    self.runner    = ClusterJobRunner()
    self.submitter = ClusterJobSubmitter()

    self.memsize   = memsize   # expressed in GB
    self.throttle  = throttle
    self.ssh_key   = ssh_key

    local_workdir = os.path.abspath(local_workdir)
    if not os.path.exists(local_workdir):
      os.mkdir(local_workdir)
    self.local_workdir = local_workdir

  def submit_command(self, cmd, *args, **kwargs):
    '''
    Submit a command to be run via bsub on the cluster. Returns the ID
    of the launched job. To wait on the completion of the submitted
    job, see the wait_on_cluster method.
    '''
    return self.submitter.submit_command(cmd, *args, **kwargs)

  def run_command(self, cmd, *args, **kwargs):
    '''
    Run a command directly on the cluster head node, waiting for the
    result to be returned. Returns a file descriptor containing the
    stdout of the job. This method is typically used for commands
    which should complete almost immediately.
    '''
    return self.runner.run_command(cmd, *args, **kwargs)

  def cluster_file_exists(self, file):
    '''
    Test whether a file exists in the configured cluster working directory.
    '''
    cmd = ("if [ -e %s ]; then echo yes; else echo no; fi" % file)
    LOGGER.debug(cmd)

    # This runs the test command in the cluster working directory.
    with self.runner.run_command(cmd) as ofh:
      first_line = ofh.readline()
      first_line = first_line.rstrip('\n')
      
    return first_line == 'yes'

  def cluster_jobs_count(self):
    '''
    Return a count of the jobs currently running on the cluster for
    the configured user.
    '''
    cmd = ("bjobs -u %s" % self.config.clusteruser)
    LOGGER.debug(cmd)
    count = 0

    with self.runner.run_command(cmd) as ofh:
      for line in ofh:
        count += 1

    # Account for the extra header line.
    return count - 1
  
  def return_file_to_localhost(self, clusterout, outfile, execute=True, donefile=False):
    '''
    If execute is False, returns a command string that can be used to
    transfer a cluster output files back to our local working
    directory. If execute is True, the command will also be run on the
    cluster.
    '''
    myhost = getfqdn()
    myuser = getuser()
    sshcmd = "scp"

    # Transferring the files back to localhost requires an appropriate
    # passwordless ssh key to be given access on our localhost. The
    # alternative is some horrendous pexpect hack which is only a
    # little more secure (see: sshSangerTunnel.py).
    if self.ssh_key is not None:
      sshcmd += " -i %s" % self.ssh_key

    # Note that we need quoting of e.g. file paths containing
    # spaces. Also, the initial './' allows filenames to contain
    # colons.
    if not os.path.isabs(clusterout):
      clusterout = './%s' % (clusterout,)
    sshcmd += (r' %s %s@%s:\"' % (clusterout, myuser, myhost)
               + bash_quote(bash_quote(self.local_workdir)) + r'/%s\"' % outfile)

    if donefile:
      sshcmd += " && ssh"
      if self.ssh_key is not None:
        sshcmd += " -i %s" % self.ssh_key
      sshcmd += (r' %s@%s touch ' % (myuser, myhost)
                 + bash_quote(bash_quote(self.local_workdir)) + r'/%s.done' % outfile)

    if execute is True:
      # This *should* die on failure.
      self.runner.run_command(sshcmd)

    return sshcmd

  def wait_on_cluster(self, jobs, cleanup_cmd=None):
    '''
    Wait for the alignment jobs running on the cluster to contact a
    designated socket file location to indicate that the jobs have
    finished.
    '''
    # Set up a job to notify localhost that the cluster is finished.
    with NamedTemporaryFile() as sobj:
      socketfile = sobj.name

    # The nc utility is pretty commonly installed; if it is not, this
    # will not work.
    LOGGER.info("Submitting monitor job to the cluster.")
    cmd = "ssh"
    if self.ssh_key is not None:
      cmd += " -i %s" % self.ssh_key
    cmd += (" %s@%s 'echo OK | nc -U %s'"
            % (getuser(), getfqdn(), socketfile))
    monjob = self.submitter.submit_command(cmd, depend_jobs=jobs,
                                         auto_requeue=False)

    # Optional clean-up job, typically used to delete temporary files.
    if cleanup_cmd is not None:
      LOGGER.info("Submitting clean-up job to the cluster.")
      self.submitter.submit_command(cleanup_cmd,
                                    depend_jobs=[monjob],
                                    auto_requeue=False)

    # Set up a socket server and wait for the cluster to get back to us.
    LOGGER.info("Waiting on a reply from the cluster...")
    sock = socket(AF_UNIX, SOCK_STREAM)
    sock.bind(socketfile)
    sock.listen(1)
    (conn, _addr) = sock.accept()

    message = ''
    while 1:
      data = conn.recv(1024)
      if not data:
        break
      message += data
    conn.close()
    os.unlink(socketfile)

    LOGGER.info("Cluster reply received: %s", message)

    return

