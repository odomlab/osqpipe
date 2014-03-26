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
from tempfile import gettempdir

from utilities import call_subprocess, bash_quote, \
    is_zipped, set_file_permissions
from config import Config

from setup_logs import configure_logging
LOGGER = configure_logging('bwa_runner')

##############################################################################

def make_bam_name(fqname):

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
            auto_requeue=False, depend_jobs=None, *args, **kwargs):
    # Pass the PYTHONPATH to the cluster process. This allows us to
    # isolate e.g. a testing instance of the code from production.
    # Note that we can't do this as easily for PATH itself because
    # bsub itself is in a custom location on the cluster.
    python_path = os.environ['PYTHONPATH']
    if not python_path:
      python_path = ''

    cmd = super(BsubCommand, self).build(cmd, *args, **kwargs)

    # Note that if this gets stuck in an infinite loop you will need
    # to use "bkill -r" to kill the job on LSF. N.B. exit code 139 is
    # a core dump. But so are several other exit codes; add 128 to all
    # the unix signals which result in a dump ("man 7 signal") for a
    # full listing.
    qval = '-Q "all ~0"' if auto_requeue else ''

    bsubcmd = (("PYTHONPATH=%s bsub -R 'rusage[mem=%d]' -r"
           + " -o %s/%%J.stdout -e %s/%%J.stderr %s")
           % (python_path,
              mem,
              self.conf.clusterstdoutdir,
              self.conf.clusterstdoutdir,
              qval))

    if queue is not None:
      bsubcmd += ' -q %s' % queue

    if jobname is not None:
      bsubcmd += ' -J %s' % jobname

    if depend_jobs is not None:
      depend = "&&".join([ "ended(%d)" % (x,) for x in depend_jobs ])
      bsubcmd += ' -w "%s"' % depend

    # To group things in a pipe (allowing e.g. use of '&&'), we use a
    # subshell. Note that we quote the sh -c string once, and
    # everything within that string twice. Commands can be of the following form:
    #
    # "ssh blah@blah.org 'echo stuff | target.txt'"
    # r"ssh blah@blah.org \"echo stuff | target.txt\""
    #
    # I.e., one needs to be careful of python's rather idiosyncratic
    # string quoting rules, and use the r"" form where necessary.
    bsubcmd += r' sh -c \"(%s)\"' % re.sub(r'"', r'\\\"', cmd)

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
  jr.submit_command(cmd, mem=10000, queue='dolab')

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

  def run_command(self, cmd, tmpdir=None, path=None, *args, **kwargs):

    cmd = self.command_builder.build(cmd, *args, **kwargs)

    if path is None:
      path = self.config.hostpath

    if tmpdir is None:
      tmpdir = gettempdir()

    LOGGER.info(cmd)
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

class RemoteJobRunner(JobRunner):
  '''
  Abstract base class holding some common methods used by classes
  controlling alignment job submission to the cluster and to other
  computing resources.
  '''
  remote_host = None
  remote_user = None
  remote_wdir = None

  def __init__(self, *args, **kwargs):

    # A little programming-by-contract, as it were.
    if not all( x in self.__dict__.keys()
                for x in ('remote_host', 'remote_user', 'remote_wdir')):
      raise StandardError("Remote host information not provided.")
    super(RemoteJobRunner, self).__init__(*args, **kwargs)

  def run_command(self, cmd, wdir=None, path=None, *args, **kwargs):
    '''
    Method used to run a command *directly* on the remote host. No
    attempt will be made to use any kind of queuing or backgrounding
    mechanism.

    The command call is wrapped in an ssh connection. This command will
    also automatically change to the configured remote working
    directory before executing the command.
    '''
    cmd = self.command_builder.build(cmd, *args, **kwargs)
    
    if wdir is None:
      wdir = self.remote_wdir

    if path is None:
      pathdef = ''
    else:
      if type(path) is list:
        path = ":".join(path)
      pathdef = "PATH=%s" % path

    cmd = ("ssh %s@%s \"cd %s && %s %s\""
           % (self.remote_user,
              self.remote_host,
              wdir,
              pathdef,
              cmd))
    LOGGER.info(cmd)
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

      destfile = os.path.join(self.remote_wdir, destfn)
      destfile = bash_quote(destfile)

      cmd = " ".join(('scp', '-p', '-q', bash_quote(fromfn),
                      "%s@%s:%s" % (self.remote_user,
                                    self.remote_host,
                                    quote(destfile)))) # FIXME sshfs?

      LOGGER.info(cmd)
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
    destfile = bash_quote(destfile)

    # Assumes that gzip is in the executable path on the remote server.
    cmd = " ".join(('gzip -f -d', quote(destfile)))
    self.run_command(cmd)

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
    self.remote_user = self.conf.clusteruser
    self.remote_wdir = self.conf.clusterworkdir if remote_wdir is None else remote_wdir

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
    for line in pout:
      matchobj = jobid_pattern.search(line)
      if matchobj:
        return int(matchobj.group(1))

    raise ValueError("Unable to parse bsub output for job ID.")

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
    self.remote_user = self.conf.althostuser
    self.remote_wdir = self.conf.althostworkdir

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
  __slots__ = ('finaldir', 'genome', 'job', 'conf')

  job = None
  
  def __init__(self, genome, finaldir='.', *args, **kwargs):

    # A little programming-by-contract, as it were.
    if not all( x in self.__dict__.keys()
                for x in ('job')):
      raise StandardError("JobRunner instance not set.")

    self.conf = Config()

    # Support relative paths as input.
    self.finaldir = os.path.realpath(finaldir)

    # Check if genome exists.
    LOGGER.info("Checking if specified genome file exists.")
    cmd = ("if [ -f %s ]; then echo yes; else echo no; fi" % genome)
    LOGGER.debug(cmd)

    if not self.job.test_mode:
      cmdstdoutfile = self.job.run_command(cmd)
      first_line = cmdstdoutfile.readline()
      first_line = first_line.rstrip('\n')
      if first_line != 'yes':
        raise ValueError("Genome %s unacessible or missing." % genome)

    self.genome = genome

  @classmethod
  def genome_path(cls, genome, indexdir, genomedir):
    '''
    Returns the expected path to the fasta file for a given genome.
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

    # Next, submit the actual jobs on the actual cluster.
    if is_paired:
      LOGGER.debug("Running bwa on paired-end sequencing input.")
      fnlist = " ".join([ quote(x) for x in destnames ])
      ## FIXME think about ways this could be improved.
      ## In the submitted command:
      ##   --rcp       is where cs_runBwaWithSplit_Merge.py eventually copies
      ##                 the reassembled bam file (via scp).
      cmd = ("%s --loglevel %d %s %s --rcp %s:%s %s %s"
             % ('cs_runBwaWithSplit.py',
                LOGGER.getEffectiveLevel(),
                cleanupflag,
                noccflag,
                self.conf.datahost,
                self.finaldir,
                self.genome,
                fnlist))

    else:
      LOGGER.debug("Running bwa on single-end sequencing input.")
      fnlist = quote(destnames[0])
      cmd = ("%s --loglevel %d %s %s --rcp %s:%s %s %s"
             % ('cs_runBwaWithSplit.py',
                LOGGER.getEffectiveLevel(),
                cleanupflag,
                noccflag,
                self.conf.datahost,
                self.finaldir,
                self.genome,
                fnlist))

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
    destnames = self.job.transfer_data(filenames, destnames)

    outfnbase = make_bam_name(destnames[0])
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

    # This is invariant PE vs. SE.

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
                                 ssh_path=conf.althostpath)

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

class BwaRunner(object):

  '''Parent class handling various functions required by scripts
  running BWA across multiple parallel alignments (and merging their
  output) on the cluster.'''

  __slots__ = ('conf', 'bwa_prog', 'samtools_prog',
               'merge_prog', 'logfile', 'debug')

  def __init__(self, debug=True):

    self.conf = Config()

    # These are now identified by passing in self.conf.clusterpath to
    # the remote command.
    self.bwa_prog      = 'bwa'
    self.samtools_prog = 'samtools'
    self.merge_prog    = 'cs_runBwaWithSplit_Merge.py'
    self.logfile       = self.conf.splitbwarunlog
    self.debug         = debug

  def _configure_logging(self, name, logger=LOGGER):
    """Configures the logs to be saved in self.logfile"""

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

  @staticmethod
  def set_file_permissions(group, path):
    """Sets group for path"""
    gid = grp.getgrnam(group).gr_gid
    os.chown(path, -1, gid)
    os.chmod(path,
             stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IWGRP|stat.S_IROTH)

##############################################################################

class SplitBwaRunner(BwaRunner):

  '''Class used to launch the initial file splitting and bwa
  alignments. This class also submits a job dependent on the outputs
  of those alignments, which in turn merges the outputs to generate
  the final bam file.'''

  # N.B. we allow creation of a __dict__ in this class since __slots__
  # is already used in the superclass.

  def __init__(self, cleanup=False, loglevel=logging.WARNING, reads=1000000,
                group=None, nocc=None, *args, **kwargs):
    super(SplitBwaRunner, self).__init__(*args, **kwargs)
    self.cleanup = cleanup
    self.reads   = reads
    self.group   = group
    self._configure_logging(self.__class__.__name__, LOGGER)
    LOGGER.setLevel(loglevel)

    if nocc:
      self.nocc = "-n %s" % (nocc,)
    else:
      self.nocc = ''

    self.bsub=JobRunner(command_builder=BsubCommand())
    # Each of these is called from within a python process which has
    # itself had PATH set appropriately.
    self.commands = {
      'SPLIT'       : "split -l %s %s %s", # split -l size file.fq prefix
      'BWA_SE'      : "%s aln %s %s | %s samse %s %s - %s"
                                + " | %s view -b -S -u - > %s",
      'BWA_PE1'     : "%s aln %s %s > %s",
      'BWA_PE2'     : "%s sampe %s %s %s %s %s %s | %s view -b -S -u - > %s",
      'MERGE'       : "%s --loglevel %d %s %s %s %s",
      'MERGE_RCP'   : "%s --loglevel %d %s %s --rcp %s %s %s",
      }

  def split_fq(self, fastq_fn):
    """Splits fastq file to 1M read per file using linux command line split"""

    LOGGER.debug("splitting fq file %s to %s per file ", fastq_fn, self.reads)

    fastq_fn_suffix = fastq_fn + '-'
    cmd = self.commands['SPLIT'] % (self.reads*4, quote(fastq_fn),
                                    quote(fastq_fn_suffix))
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

  def _submit_lsfjob(self, command, jobname, depend=None):
    """ Executes command in LSF cluster """

    jobid = self.bsub.submit_command(command, jobname=jobname,
                                     depend_jobs=depend, mem=10000,
                                     path=self.conf.clusterpath,
                                     tmpdir=self.conf.clusterworkdir,
                                     queue=self.conf.clusterqueue)
    return '' if jobid is None else jobid

  def run_bwas(self, genome, paired, fq_files, fq_files2):
    """Submits bwa alignment jobs for list of fq files to LSF cluster"""

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
        cmd = self.commands['BWA_PE1'] % (self.bwa_prog, genome,
                             bash_quote(fqname),
                             bash_quote(sai_file1))
        cmd2 = self.commands['BWA_PE1'] % (self.bwa_prog, genome,
                              bash_quote(fq_files2[current]),
                              bash_quote(sai_file2))
        cmd3 = self.commands['BWA_PE2'] % (self.bwa_prog, self.nocc, genome,
                              bash_quote(sai_file1),
                              bash_quote(sai_file2),
                              bash_quote(fqname),
                              bash_quote(fq_files2[current]),
                              self.samtools_prog, out)

        LOGGER.info("starting bwa step1 on '%s'", fqname)
        jobid_sai1 = self._submit_lsfjob(cmd, jobname1)
        LOGGER.debug("got job id '%s'", jobid_sai1)
        LOGGER.info("starting bwa step1 on '%s'", fq_files2[current])
        jobid_sai2 = self._submit_lsfjob(cmd2, jobname2)
        LOGGER.debug("got job id '%s'", jobid_sai2)

        if jobid_sai2 and jobid_sai2:
          LOGGER.info("preparing bwa step2 on '%s'", fqname)
          jobid_bam = self._submit_lsfjob(cmd3, jobname_bam, (jobid_sai1, jobid_sai2))
          LOGGER.debug("got job id '%s'", jobid_sai2)
          job_ids.append(jobid_bam)
        else:
          LOGGER.error("bjob submission for bwa step1 for '%s' or '%s' failed!",
                       fqname, fq_files2[current])
      else:
        cmd = self.commands['BWA_SE'] % (self.bwa_prog, genome,
                            bash_quote(fqname),
                            self.bwa_prog, self.nocc,
                            genome,
                            bash_quote(fqname),
                            self.samtools_prog, out)
        LOGGER.info("starting bwa on '%s'", fqname)
        LOGGER.debug(cmd)
        jobid_bam = self._submit_lsfjob(cmd, jobname_bam)
        LOGGER.debug("got job id '%s'", jobid_bam)
        job_ids.append(jobid_bam)
      current += 1

    return (job_ids, out_names)

  def queue_merge(self, bam_files, depend, bam_fn, rcp_target):
    """Submits samtools job for merging list of bam files to LSF cluster"""
    input_files = " ".join(bam_files) # singly-bash-quoted
    LOGGER.debug("Entering queue_merge with input_files=%s", input_files)
    cleanup = ""
    group = ""
    cmd = ""
    jobname = bam_files[0].split("_")[0] + "bam"

    if self.cleanup:
      cleanup = "--cleanup"

    if self.group:
      group = "--group %s" % (self.group,)

    if rcp_target:
      cmd = self.commands['MERGE_RCP'] % (self.merge_prog,
                                          LOGGER.getEffectiveLevel(),
                                          cleanup, group, rcp_target,
                                          bash_quote(bam_fn), input_files)
    else:
      cmd = self.commands['MERGE'] % (self.merge_prog,
                                      LOGGER.getEffectiveLevel(),
                                      cleanup, group,
                                      bash_quote(bam_fn), input_files)

    LOGGER.info("Preparing bwa merge on '%s'", input_files)
    LOGGER.debug(cmd)
    jobid = self._submit_lsfjob(cmd, jobname, depend)
    LOGGER.debug("got job id '%s'", jobid)

  def run(self, files, genome, rcp_target=None):

    '''Main entry point for the class.'''

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
    bam_fn = make_bam_name(files[0])

    self.queue_merge(bam_files, job_ids, bam_fn, rcp_target)

##############################################################################

class MergeBwaRunner(BwaRunner):

  '''Class used to merge a set of bam files into a single output bam
  file.'''

  def __init__(self, cleanup=False, loglevel=logging.WARNING,
               group=None, *args, **kwargs):
    super(MergeBwaRunner, self).__init__(*args, **kwargs)
    self.cleanup = cleanup
    self.group   = group
    self._configure_logging(self.__class__.__name__, LOGGER)
    self.commands = {
      'MERGE' : "%s merge - %s | %s sort - %s",
      'COPY'  : "scp -p -q %s %s",
      'DONE'  : "ssh %s touch %s/%s.done",
      }
    LOGGER.setLevel(loglevel)

  @staticmethod
  def input_files_exist(input_fns):
    """Checks if input files exist"""
    okay = True
    for fname in input_fns:
      if not os.path.exists(fname):
        LOGGER.error("missing expected file '%s', cannot continue.", fname)
        okay = False
      elif not os.path.isfile(fname):
        LOGGER.error("file '%s' is not a regular file, cannot continue.", fname)
        okay = False
    return okay

  def merge_files(self, output_fn, input_fns):
    """Merges list of bam files"""
    output_fnfull = output_fn + ".bam"
    if len(input_fns) == 1:
      LOGGER.warn("renaming file: %s", input_fns[0])
      move(input_fns[0], output_fnfull)
    else:
      cmd = (self.commands['MERGE'] % (self.samtools_prog,
                          " ".join([ bash_quote(x) for x in input_fns]),
                          self.samtools_prog,
                          bash_quote(output_fn)))
      LOGGER.debug(cmd)
      pout = call_subprocess(cmd, shell=True,
                            tmpdir=self.conf.clusterworkdir,
                            path=self.conf.clusterpath)
      for line in pout:
        LOGGER.warn("SAMTOOLS: %s", line[:-1])
    if not os.path.isfile(output_fnfull):
      LOGGER.error("expected output file '%s' cannot be found.", output_fnfull)
      sys.exit("File access error.")
    if self.group:
      self.set_file_permissions(self.group, output_fn)
    if self.cleanup:
      for fname in input_fns:
        # remove bam file
        if len(input_fns) > 1:
          LOGGER.info("Unlinking bam file '%s'", fname)
          os.unlink(fname)
        # strip bam to fastq and remove
        fq1 = os.path.splitext(fname)[0]
        LOGGER.info("Unlinking fastq file '%s'", fq1)
        os.unlink(fq1)
        # strip fastq to base and check if paired end
        (fq1base, fq1ext) = os.path.splitext(fq1)
        base_pattern = re.compile(r'^(.*\d\d)(p\d)$')
        matchobj = base_pattern.match(fq1base)
        # if paired end pattern found, remove pe fastq and sai files
        if matchobj:
          sext = matchobj.group(2)
          if sext == 'p1':
            pext = 'p2'
          elif sext == 'p2':
            pext = 'p1'
          else:
            raise ValueError("Unexpected paired-end designation: %s" % sext)
          fq2 = matchobj.group(1) + pext + fq1ext
          sai1 = fq1 + ".sai"
          sai2 = fq2 + ".sai"
          LOGGER.info("Unlinking fastq file '%s'", fq2)
          os.unlink(fq2)
          LOGGER.info("Unlinking sai file '%s'", sai1)
          os.unlink(sai1)
          LOGGER.info("Unlinking sai file '%s'", sai2)
          os.unlink(sai2)

  def copy_result(self, target, fname):
    """Copies file to target location"""
    qname = bash_quote(fname)
    cmd = self.commands['COPY'] % (qname, target)
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
      cmd = self.commands['DONE'] % (flds[0], flds[1], bash_quote(fn_base))
      LOGGER.debug(cmd)
      call_subprocess(cmd, shell=True,
                     tmpdir=self.conf.clusterworkdir,
                     path=self.conf.clusterpath)
    if self.cleanup:
      os.unlink(fname)
    return

  def run(self, input_fns, output_fn, rcp_target=None):
    '''Main entry point for the class.'''
    LOGGER.info("merging '%s' into '%s'", ", ".join(input_fns), output_fn)
    self.merge_files(output_fn, input_fns)
    LOGGER.info("merged '%s' into '%s'", ", ".join(input_fns), output_fn)
    if rcp_target:
      self.copy_result(rcp_target, output_fn + ".bam")
      LOGGER.info("copied '%s' to '%s'", output_fn, rcp_target)

