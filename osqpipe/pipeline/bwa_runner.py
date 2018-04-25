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
Parent class for the code running BWA on the cluster. Here we
handle some core external dependencies, logging and group
permissions.
'''

import sys
import os
import re
import stat
import grp
from pipes import quote
from distutils import spawn
from socket import getfqdn, socket, AF_UNIX, SOCK_STREAM, gethostname
from getpass import getuser
from tempfile import NamedTemporaryFile

from osqutil.utilities import bash_quote, sanitize_samplename, \
  BamPostProcessor
from osqutil.progsum import ProgramSummary
from osqutil.cluster import make_bam_name_without_extension, \
  ClusterJobRunner, ClusterJobSubmitter, DesktopJobSubmitter
from osqutil.config import Config

from osqutil.setup_logs import configure_logging
LOGGER = configure_logging('bwa_runner')

##############################################################################

def genome_fasta_path(genome, genomedir, indexdir=None):
  '''
  Returns the expected path to the fasta file for a given genome
  index. If a specific fasta index directory is specified
  (e.g. bwa-0.6.1), the path points to a fasta file in that
  subdirectory.
  '''
  sciname = genome.species.scientific_name
  sciname = sciname.replace(" ", "_")
  sciname = sciname.lower()
  if indexdir is not None:
    gpath = os.path.join(genomedir, sciname,
                         genome.code, indexdir, "%s.fa" % genome.code)
  else:
    gpath = os.path.join(genomedir, sciname,
                         genome.code, "%s.fa" % genome.code)

  return gpath

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

class AlignmentJobRunner(object):
  '''
  An abstract class from which all alignment job runners inherit. All
  subclasses must instantiate a RemoteJobRunner class in the 'job' slot
  prior to calling the base __init__ method.
  '''

  # ML: Note that slot 'aligner' has been introduced only for the sake of 'star' aligner
  # as the method for checking presence for reference genome index is different.
  # We may want to find better way of dealing with this.
  __slots__ = ('finaldir', 'genome', 'job', 'conf', 'samplename', 'aligner')

  job = None
  
  def __init__(self, genome, finaldir='.', samplename=None, aligner=None, *args, **kwargs):

    # A little programming-by-contract, as it were.
#    if not all( hasattr(self, x) for x in ('job')):
#      raise StandardError("JobRunner instance not set.")

    self.conf = Config()

    # Support relative paths as input.
    self.finaldir = os.path.realpath(finaldir)

    # Check if genome exists.
    LOGGER.info("Checking if specified genome file exists.")
    cmd = None
    if aligner is not None and aligner=='star':
      cmd = ("if [ -d %s ]; then echo yes; else echo no; fi" % genome)
    else:
      cmd = ("if [ -f %s ]; then echo yes; else echo no; fi" % genome)
    LOGGER.debug(cmd)

    if not self.job.test_mode:
      runjob = ClusterJobRunner(test_mode=self.job.test_mode)
      cmdstdoutfile = runjob.run_command(cmd)
      first_line = cmdstdoutfile.readline()
      first_line = first_line.rstrip('\n')
      if first_line != 'yes':
        raise ValueError("Genome %s inaccessible or missing." % genome)

    self.genome = genome
    self.samplename = sanitize_samplename(samplename)

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
             nocc=None, bwa_algorithm='aln', fileshost=None, nosplit=False, rcp=None, lcp=None, *args, **kwargs):

    '''Actually submit the job. The optional destnames argument can be
    used to name files on the cluster differently to the source. This
    is occasionally useful.'''

    assert(bwa_algorithm in ('aln', 'mem'))
    paired_sanity_check(filenames, is_paired)

    # by lukk01:
    # NB! Copying files to cluster is not any more necessary as long s the hostflag = '--fileshost %s' is uncommented below.
    # However, this would be a pull rather than push and we should then pull from the archive for process_file.py
    # It would requre, though, re-writing of data processing and alignment orders in process_file.py. Just a throught.
    #
    
    # First, copy the files across and uncompress on the server.
    LOGGER.info("Copying files to the cluster.")
    destnames = self.job.transfer_data(filenames, destnames)

    # Next, create flag for cleanup
    cleanupflag = '--cleanup' if cleanup else ''

    # Next, create flag for number of non-unique reads to keep in samse/sampe
    noccflag = ('--n_occ %s' % (nocc,)) if nocc else ''

    # Sample names containing spaces are bad on the command line,
    # and potentially problematic in bam read groups.
    sampleflag = '--sample %s' % self.samplename if self.samplename else ''

    # Whether to run bwa mem or aln.
    algoflag = '--algorithm %s' % bwa_algorithm

    # Deal with default values for fileshost and rcp/lcp. I.e. figure out if files are located in cluster and results would need to be copied somewhere or not.
    cpflag = '' 
    hostflag = ''
    filehost = gethostname()
    if filehost != self.conf.cluster:
      # hostflag  = '--fileshost %s' % filehost
      cpflag = '--rcp %s:%s' % (self.conf.datahost, self.finaldir)
    else:
      # the files are already in host. Override cleanup to prevent source files to be deleted.
      LOGGER.info("Input files are local. Overriding --cleanup to prevent files being deleted.")
      cleanupflag = ''
      cpflag = '--lcp %s' % self.finaldir

    # If fileshost has been specified, override default
    if fileshost is not None:
      hostflag = '--fileshost %s' % fileshost
    # If rcp has been specified, override default
    if rcp is not None:
      cpflag = '--rcp %s' % rcp
    # If lcp has been specified, override default
    if lcp is not None:
      cpflag = '--lcp %s' % lcp
    # If nosplit has been set, forward the value
    splitflag = ''
    if nosplit is not None:
      splitflag = '--no-split'

    # This now searches directly on the cluster.
    progpath = self.job.find_remote_executable('cs_runBwaWithSplit.py',
                                               path=self.conf.clusterpath)

    if progpath is None:
      raise StandardError("cs_runBwaWithSplit.py not found on clusterpath. Possible misconfiguration?")

    # Next, submit the actual jobs on the actual cluster.
    if is_paired:
      LOGGER.debug("Running bwa on paired-end sequencing input.")
      fnlist = " ".join([ quote(x) for x in filenames ])
      # fnlist = " ".join([ quote(x) for x in destnames ])
      ## FIXME think about ways this could be improved.
      ## In the submitted command:
      ##   --rcp       is where cs_runBwaWithSplit_Merge.py eventually copies
      ##                 the reassembled bam file (via scp).
      cmd = ("python %s --loglevel %d %s %s %s %s %s %s %s %s %s"
             % (progpath,
                LOGGER.getEffectiveLevel(),
                cleanupflag,
                hostflag,
                noccflag,
                cpflag,
                splitflag,
                sampleflag,
                algoflag,
                self.genome,
                fnlist))

    else:
      LOGGER.debug("Running bwa on single-end sequencing input.")
      fnlist = quote(filenames[0])
      # fnlist = quote(destnames[0])
      cmd = ("python %s --loglevel %d %s %s %s %s %s %s %s %s %s"
             % (progpath,
                LOGGER.getEffectiveLevel(),
                cleanupflag,
                hostflag,
                noccflag,
                cpflag,
                splitflag,
                sampleflag,
                algoflag,
                self.genome,
                fnlist))

    LOGGER.info("Submitting bwa job to cluster.")
    self.job.submit_command(cmd, *args, **kwargs)

  @classmethod
  def build_genome_index_path(cls, genome, *args, **kwargs):

    # Import here rather than main file as otherwise cluster operations fail.
    from ..models import Program

    conf = Config()

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

    gpath = genome_fasta_path(genome, indexdir=indexdir, genomedir=conf.clustergenomedir)

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

    # First, copy the files across and uncompress on the server. We
    # remove commas here because otherwise tophat is a little too keen
    # to split on them (quoting doesn't work).
    LOGGER.info("Copying files to the cluster.")
    destnames = [ re.sub(',+', '_', os.path.basename(fname)) for fname in filenames ]
    destnames = self.job.transfer_data(filenames, destnames)

    # Next, create flag for cleanup
    if cleanup:
      cleanupflag = '--cleanup'
    else:
      cleanupflag = ''

    if self.samplename:
      sampleflag = '--sample %s' % self.samplename
    else:
      sampleflag = ''

    # This now searches directly on the cluster.
    progpath = self.job.find_remote_executable('cs_runTophatWithSplit.py',
                                               path=self.conf.clusterpath)

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

    # Import here rather than main file as otherwise cluster operations fail.
    from ..models import Program

    conf = Config()

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
    gpath = genome_fasta_path(genome, indexdir=indexdir, genomedir=conf.clustergenomedir)

    return gpath

##############################################################################

class StarClusterJobSubmitter(AlignmentJobRunner):
  '''
  Class representing the submission of a STAR job to the
  cluster. This class works similarly to BwaClusterJobSubmitter.
  '''
  def __init__(self, *args, **kwargs):
    self.job = ClusterJobSubmitter(*args, **kwargs)
    super(StarClusterJobSubmitter, self).__init__(*args, **kwargs)

  def submit(self, filenames, is_paired=False, destnames=None, cleanup=True,
              *args, **kwargs):
    '''
    Actually submit the job. The optional destnames argument can be
    used to name files on the cluster differently to the source. This
    is occasionally useful.
    '''
    paired_sanity_check(filenames, is_paired)

    # First, copy the files across and uncompress on the server. We
    # remove commas here because otherwise tophat is a little too keen
    # to split on them (quoting doesn't work).
    LOGGER.info("Copying files to the cluster.")
    destnames = [ re.sub(',+', '_', os.path.basename(fname)) for fname in filenames ]
    destnames = self.job.transfer_data(filenames, destnames)

    # Next, create flag for cleanup
    if cleanup:
      cleanupflag = '--cleanup'
    else:
      cleanupflag = ''

    if self.samplename:
      sampleflag = '--sample %s' % self.samplename
    else:
      sampleflag = ''

    # This now searches directly on the cluster.
    progpath = self.job.find_remote_executable('cs_runStarWithSplit.py',
                                               path=self.conf.clusterpath)

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

    LOGGER.info("Submitting STAR job to cluster.")
    self.job.submit_command(cmd, *args, **kwargs)

  @classmethod
  def build_genome_index_path(cls, genome, *args, **kwargs):

    # Import here rather than main file as otherwise cluster operations fail.
    from ..models import Program

    conf = Config()

    # Get information about default aligner, check that the program is
    # in path and try to predict its version.
    alignerinfo = ProgramSummary('STAR',
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
      indexdir = "%s_%s" % ('STAR', alignerinfo.version)

    except Program.DoesNotExist, _err:
      sys.exit(("""Aligner "%s" version "%s" found at path "%s" """
               % (alignerinfo.program, alignerinfo.version, alignerinfo.path))
               + "not recorded as current in repository! Quitting.")

    # Build path to STAR genome dir. Note that STAR takes dir name only without indexdir.fa suffix in the end.
    gpath = genome_fasta_path(genome, indexdir=indexdir, genomedir=conf.clustergenomedir)
    # A bit of an ugly hack here: Remove indexdir.fa suffix from gpath created by genome_fasta_path
    gpath = os.path.split(gpath)[0]
    
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
    self.num_threads   = num_threads
    self.tempfiles     = []

  def _run_pairedend_bwa_aln(self, destnames, outfnbase, noccflag=''):
    
    LOGGER.debug("Running bwa aln on paired-end sequencing input.")
    num_threads = max(1, int(self.num_threads))
    fnlist = " ".join([ quote(x) for x in destnames ])

    # We have to run bwa aln on both inputs separately, then combine
    # them with bwa sampe.
    fq_0 = quote(destnames[0])
    fq_1 = quote(destnames[1])

    sai_0 = "%s.sai" % fq_0
    sai_1 = "%s.sai" % fq_1
    self.tempfiles.extend([ sai_0, sai_1 ])

    # FIXME note that the PATH is not set correctly for pipe sinks
    # in this command (e.g. samtools). Either we need to wrap this
    # somehow (see SplitBwaRunner) or we need to drop PATH in after
    # each '|'.
    cmd  = ("bwa aln -t %d %s %s > %s"
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
                                        + " - | samtools sort -@ %d -o %s.bam -")
            % (noccflag,
               self.genome,
               sai_0,
               sai_1,
               fq_0,
               fq_1,
               num_threads,
               num_threads,
               outfnbase))
    return cmd

  def _run_singleend_bwa_aln(self, destnames, outfnbase, noccflag):

    LOGGER.debug("Running bwa aln on single-end sequencing input.")
    num_threads = max(1, int(self.num_threads))
    fnlist = quote(destnames[0])
    cmd  = (("bwa aln -t %d %s %s | bwa samse %s %s - %s"
                              + " | samtools view -b -S -u -@ %d -"
                              + " | samtools sort -@ %d -o %s.bam -")
            % (num_threads,
               self.genome,
               fnlist,
               noccflag,
               self.genome,
               fnlist,
               num_threads,
               num_threads,
               outfnbase))

    return cmd
  
  def _run_bwa_mem(self, destnames, outfnbase, noccflag):

    nfq = len(destnames)
    if nfq == 1:
      LOGGER.debug("Running bwa mem on single-end sequencing input.")
    elif nfq == 2:
      LOGGER.debug("Running bwa mem on paired-end sequencing input.")
    else:
      raise StandardError("Incorrect number of files passed to bwa mem: %d" % nfq)
    
    num_threads = max(1, int(self.num_threads))
    sortstr = "samtools sort -@ %d -o %s.bam -"
    fnlist = " ".join([quote(fqname) for fqname in destnames ])
    cmd  = (("bwa mem -t %d %s %s | samtools view -b -S -u -@ %d -"
             + " | %s")
            % (num_threads,
               self.genome,
               fnlist,
               num_threads,
               num_threads,
               outfnbase, sortstr))
    
    return cmd
  
  def submit(self, filenames,
             is_paired=False, destnames=None, cleanup=True,
             nocc=None, bwa_algorithm='aln', *args, **kwargs):

    '''Submit a job as a background chained process on the designated
    remote server.'''

    assert(bwa_algorithm in ('aln', 'mem'))
    paired_sanity_check(filenames, is_paired)

    # First, copy the files across and uncompress on the server.
    LOGGER.info("Copying files to the alignment server.")
    destnames = self.job.transfer_data(filenames, destnames)

    outfnbase = make_bam_name_without_extension(destnames[0])
    outfnfull = outfnbase + '.bam'
    outfnfullout = outfnbase + '.bam'

    compress = True
    if self.conf.compressintermediates == "False":
      compress = False
      outfnfullout = outfnbase + '_uc.bam'
      
    self.tempfiles = destnames + [ outfnfull ]

    # Next, create flag for number of non-unique reads to keep in samse/sampe
    if nocc:

      if bwa_algorithm == 'mem':
        raise StandardError("The nocc argument is not supported by bwa mem. Try bwa aln instead.")

      if is_paired:
        noccflag = '-n %s -N %s' % (nocc, nocc)
      else:
        noccflag = '-n %s' % nocc
    else:
      noccflag = ''

    cmd = ''

    if bwa_algorithm == 'aln':

      if is_paired:
        cmd += self._run_pairedend_bwa_aln(destnames, outfnbase, noccflag)

      else:
        cmd += self._run_singleend_bwa_aln(destnames, outfnbase, noccflag)

    elif bwa_algorithm == 'mem':
      cmd += self._run_bwa_mem(destnames, outfnbase, noccflag)

    else:
      raise ValueError("BWA algorithm not recognised: %s" % bwa_algorithm)

    # This is invariant PE vs. SE. First, run our standard picard cleanup:
    postproc = BamPostProcessor(input_fn=outfnfull, output_fn=outfnfullout,
                                samplename=self.samplename,
                                tmpdir=self.conf.althostworkdir, compress=compress)
    cmd += (" && %s && rm %s"
            % (" ".join(postproc.clean_sam()), outfnfull))
    cmd += (" && %s && rm %s"
            % (" ".join(postproc.add_or_replace_read_groups()), postproc.cleaned_fn))
    cmd += (" && %s && rm %s"
            % (" ".join(postproc.fix_mate_information()), postproc.rgadded_fn))
    # if interemediate files were uncompressed bams, compress the final bam
    if not compress:
      num_threads = max(1, int(self.conf.num_threads))
      cmd += (" && samtools view -b -S -@ %d -o %s %s && rm %s" %
              (num_threads, outfnfull, outfnfullout, outfnfullout))

    # We generate verbose logging here to better monitor file transfers.
    cmd += (" && scp -v -c aes128-cbc -i %s %s %s@%s:%s"
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
      cmd += (" && rm %s" % " ".join(self.tempfiles))

    LOGGER.info("Submitting bwa job to alignment host.")
    self.job.submit_command(cmd)

  @classmethod
  def build_genome_index_path(cls, genome, *args, **kwargs):

    # Import here rather than main file as otherwise cluster operations fail.
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

    gpath = genome_fasta_path(genome, indexdir=indexdir, genomedir=conf.althostgenomedir)

    return gpath

##############################################################################
##############################################################################

class ClusterJobManager(object):
  '''
  Moderately abstract base class providing some methods and attributes
  commonly used in higher-level cluster process management classes
  (e.g. GsnapManager, LastzManager).
  '''
  __slots__ = ('namespace', 'submitter', 'runner', 'config',
               'throttle', 'memsize', 'ssh_key', 'local_workdir', 'time_limit')

  def __init__(self, namespace=None, throttle=0, memsize=20, time_limit=48,
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
    self.time_limit = time_limit
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
    sshcmd += (r' %s %s@%s:\"' % (bash_quote(clusterout), myuser, myhost)
               + bash_quote(bash_quote(self.local_workdir + r'/%s' % outfile)) + r'\"')

    if donefile:
      sshcmd += " && ssh"
      if self.ssh_key is not None:
        sshcmd += " -i %s" % self.ssh_key
      sshcmd += (r' %s@%s touch ' % (myuser, myhost)
                 + bash_quote(bash_quote(self.local_workdir + r'/%s.done' % outfile)))

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

