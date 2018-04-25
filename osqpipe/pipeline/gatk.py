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
Code used to manage the GATK cleanup parts of our HCC pipeline.
'''

import os
import re
from time import time, sleep
from subprocess import Popen, PIPE, STDOUT, CalledProcessError
from shutil import copy
from pipes import quote
from django.db import transaction

from ..models import Alnfile, Library, Alignment, MergedAlnfile, Genome, Lane
from osqutil.samtools import count_bam_reads
from osqutil.utilities import call_subprocess, checksum_file, \
    sanitize_samplename
from osqutil.config import Config
from .bwa_runner import ClusterJobManager

import xml.etree.ElementTree as ET
from tempfile import NamedTemporaryFile
from pysam import AlignmentFile
from shutil import move

from osqutil.setup_logs import configure_logging
LOGGER = configure_logging('gatk')
CONFIG = Config()

################################################################################
# A handful of utility functions.
def retrieve_readgroup_alignment(rgroup, genome=None, bamfilter=False):
  '''
  Simply returns the osqpipe Alignment object for a given read group,
  a list of which is returned by this call:

  pysam.AlignmentFile.header.get('RG', []).
  '''
  alns = Alignment.objects.filter(lane__library__code=rgroup.get('LB'),
                                  lane__facility__code=rgroup.get('CN'),
                                  lane__lanenum=rgroup.get('PU'))
  if genome is not None:
    alns = alns.filter(genome__code=genome)

  # Only include alignments with bam files. Can be useful when there are older alignments which are only present to link to earlier MergedAlignments.
  if bamfilter:
    alns = alns.filter(alnfile__filetype__code='bam')

  if alns.count() > 1:
    raise ValueError("Multiple Alignments match read group and genome parameters. Maybe filter by genome or bam file presence?")
  elif alns.count() == 0:
    raise Alignment.DoesNotExist("No Alignments found to match read group and genome parameters.")
  else:
    return alns[0]

def check_bam_readcount(bam, maln, readcountdicts=None):
  '''
  In principle, the total reads returned by the GATK pipeline should
  be the sum of the reads in the original fastq files. We check that
  here. Arguments are: bam (the GATK output bam file), maln (the
  MergedAlignment db object), readcountdicts, an (optional) tuple of
  dicts (total, mapped, munique) keyed by readgroup ID.
  '''
  expected = sum([ aln.lane.total_passedpf for aln in maln.alignments.all() ])
  if readcountdicts is None:
    numreads = count_bam_reads(bam)
  else:
    numreads = sum(readcountdicts[0].values())

  ## See how things pan out: if we have to relax this check, here
  ## would be a good place to start (i.e., raise a warning rather than
  ## an Exception).
  if numreads != expected:
    message = ("Number of reads in bam file is differs from that in "
               + "fastq file: %d (bam) vs %d (fastq)")
    raise ValueError(message % (numreads, expected))

def autocreate_alignment(rgroup, genome, readcounts):
  '''
  Given the metadata harvested from a bam file, create an Alignment
  object suitable for linking to MergedAlignment. The readcounts
  argument is a tuple or list with three entries: (total, mapped, munique).
  '''
  if genome is None:
    raise ValueError("Genome code is required for Alignment autocreation.")
  
  library = Library.objects.get(code=rgroup.get('LB'))
  lane    = Lane.objects.get(library__code=rgroup.get('LB'),
                             facility__code=rgroup.get('CN'),
                             lanenum=rgroup.get('PU'))
  genobj  = Genome.objects.get(code=genome)
  aln     = Alignment.objects.get_or_create(lane=lane,
                                            genome=genobj,
                                            total_reads=readcounts[0],
                                            mapped=readcounts[1],
                                            munique=readcounts[2])

  return aln[0]

################################################################################
# Functions to update bam read group information in as lightweight a
# manner as possible. These are dependent on samtools.
def update_bam_readgroups(bam):
  '''
  A lightweight method for updating bam read group information to
  match the current repository annotation.
  '''
  # Support both the filename or an Alnfile/MergedAlnfile object being
  # passed.
  if issubclass(type(bam), MergedAlnfile):
    _update_mergedalnfile_bam_readgroups(bam)
  elif issubclass(type(bam), Alnfile):
    _update_alnfile_bam_readgroups(bam)
  else:
    try:
      bam = Alnfile.objects.get(filename=bam)
      _update_alnfile_bam_readgroups(bam)
    except Alnfile.DoesNotExist:
      try:
        bam = MergedAlnfile.objects.get(filename=bam)
        _update_mergedalnfile_bam_readgroups(bam)
      except MergedAlnfile.DoesNotExist:
        raise ValueError("Requested bam file does not exist in database as Alnfile or MergedAlnfile: %s" % bam)

def _edit_bam_readgroup_data(bam, platform_unit=None, library=None, sample=None, center=None):
  '''
  Fairly generic internal function which makes the actual readgroup
  annotation changes. This function is deliberately agnostic about
  where the annotation comes from; it is up to the caller to make that
  decision. Note that it is assumed that the caller function is within
  a transaction; this allows us to be sure that database-derived
  annotation passed to this function will not change during the
  procedure.
  '''
  if bam.filetype.code != 'bam':
    raise ValueError("Function requires bam file, not %s (%s)"
                     % (bam.filetype.code, bam.filename))

  # First, extract the current file header.
  LOGGER.info("Reading current bam file header.")
  cmd = ('samtools', 'view', '-H', bam.repository_file_path)
  subproc = Popen(cmd, stdout=PIPE, stderr=STDOUT)
  header  = subproc.communicate()[0]
  retcode = subproc.wait()

  if retcode != 0:
    raise CalledProcessError(retcode, " ".join(cmd))

  if len(header) == 0:
    raise ValueError("The bam file has no header information; is this actually a bam file?")

  newheader = []
  for line in header.split("\n"):

    # Make the actual changes here.
    if re.match('@RG', line):
      LOGGER.info("Editing @RG header line.")
      fields = dict( field.split(':', 1) for field in line.split("\t")
                     if not re.match('^@', field) )
      if platform_unit is not None:
        fields['PU'] = platform_unit
      if library is not None:
        fields['LB'] = library
      if sample is not None:
        fields['SM'] = sample
      if center is not None:
        fields['CN'] = center
      newline = "@RG\t%s" % "\t".join([ "%s:%s" % (key, val)
                                        for (key, val) in fields.iteritems() ])
      newheader.append(newline)
    else:
      newheader.append(line)
  newheader = "\n".join(newheader)

  # Replace the old header with the edited version.
  LOGGER.info("Replacing bam file header.")
  tmpbam = "%s.reheader" % bam.repository_file_path
  move(bam.repository_file_path, tmpbam)
  cmd = 'samtools reheader - %s > %s' % (tmpbam, bam.repository_file_path)
  
  subproc = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=STDOUT, shell=True)
  stdout  = subproc.communicate(input=newheader)
  retcode = subproc.wait()

  if retcode != 0:
    raise CalledProcessError(retcode, cmd)

  LOGGER.info("Correcting bam file checksum.")
  chksum = checksum_file(bam.repository_file_path, unzip=False)
  bam.checksum = chksum
  bam.save()
  os.unlink(tmpbam)

@transaction.atomic
def _update_alnfile_bam_readgroups(bam):
  '''
  Updates read group PU, LB, SM and CN tags based on the annotation
  stored in the database linked to this Alnfile.
  '''
  bam = Alnfile.objects.get(id=bam.id)
  library = bam.alignment.lane.library
  _edit_bam_readgroup_data(bam,
                           platform_unit = bam.alignment.lane.lanenum,
                           library       = library.code,
                           sample        = sanitize_samplename(library.sample.name),
                           center        = bam.alignment.lane.facility.code)

@transaction.atomic
def _update_mergedalnfile_bam_readgroups(bam):
  '''
  Updates read group SM tag based on the annotation stored in the
  database linked to this MergedAlnfile. Note that PN, LB and CN are
  not changed, primarily because that annotation is what is used to
  link the MergedAlnfile back to its constitutent Alignments. In other
  words, if that annotation needs changing then the file has been
  mis-linked within the database, and it might be best just to remove
  the file from the repository altogether. Such cases are assumed rare
  enough that a manual fix is appropriate.
  '''
  bam = MergedAlnfile.objects.get(id=bam.id)
  samples = set([ aln.lane.library.sample.name for aln in bam.alignment.alignments.all() ])
  if len(samples) > 1:
    raise ValueError("MergedAlnfile appears to be linked to multiple samples, please fix: %s"
                     % ",".join(samples))
  _edit_bam_readgroup_data(bam,
                           sample = sanitize_samplename(list(samples)[0]))

def _make_local_bamfile_name(fname, samplename, finalprefix):
  '''
  Quick utility function so that we only construct our final bamfile
  name in one place in the code.
  '''
  return "%s%s%s" % (finalprefix, samplename, os.path.splitext(fname)[1])

################################################################################
# The main GATK handler class.
class GATKPreprocessor(ClusterJobManager):
  '''
  Class used to manage the bridge between our standard osqpipe
  pipeline and the bioinformatics core GATK preprocessing
  pipeline. This class will merge the bam files for a given set of
  libraries (which should all correspond to a single sample). It will
  then use the cluster to call picard MarkDuplicates and
  BuildBamIndex, prior to calling the GATK pipeline which is currently
  configured to run IndelRealigner and BaseRecalibrator on the input
  data, generating a bam file which is finally copied back to the
  local host.
  '''
  def __init__(self, with_markduplicates=True, with_gatkpipe=True,
               *args, **kwargs):

    super(GATKPreprocessor, self).__init__(*args, **kwargs)

    # These flags control whether this object runs MarkDuplicates
    # and/or the GATK preprocessing pipeline (IndelRealigner/BQSR).
    self.markduplicates  = with_markduplicates
    self.gatkpipe        = with_gatkpipe

    if self.gatkpipe:
      self.outprefix = 'IR_BQSR_'
    elif self.markduplicates:
      self.outprefix = 'MD_'
    else:
      self.outprefix = ''
      
  def wait_on_cluster(self, finalbam, localbam, *args, **kwargs):
    '''
    Overridden wait_on_cluster which both waits on the cluster and
    cleans up the input file if the outputs look reasonable.
    '''
    super(GATKPreprocessor, self).wait_on_cluster(*args, **kwargs)
    finaldone = "%s.done" % finalbam
    if not (os.path.exists(finalbam) and os.path.exists(finaldone)):
      raise StandardError("Expected output file %s not found, or %s not created."
                          % (finalbam, finaldone))

    # Delete local files (only if we're waiting on cluster
    # though). Note that if the cluster jobs fail, the local merged
    # bam file should not be deleted.
    LOGGER.info("Deleting local merged bam file")
    os.unlink(localbam)

  def cluster_filename(self, fname):
    '''
    Given a local filename, return a file path suitable for use on the
    cluster. This assumes that the API user is responsible for
    ensuring that the same file is not used twice within the same
    process.
    '''
    fname  = "%d_%s" % (os.getpid(), fname)
    clpath = os.path.join(CONFIG.gatk_cluster_input, fname)
    return clpath

  def samtools_merge_bams(self, bams, output_fn):
    '''
    Use samtools to merge a set of bams locally.
    '''
    if len(bams) == 0:
      raise ValueError("Zero input bam files for merging.")
    if len(bams) == 1:
      LOGGER.warning("Only one bam file supplied; copying, rather than merging.")
      copy(bams[0], output_fn)
      return
  
    # We assume our input bam files are appropriately sorted (which,
    # coming from the repository, they should be).
    cmd = ['samtools', 'merge', output_fn] + bams

    LOGGER.info("Using samtools to merge bam files: %s",
                ", ".join([ os.path.basename(bam) for bam in bams ]))

    call_subprocess(cmd, path=CONFIG.hostpath)

    if not os.path.exists(output_fn):
      raise StandardError("Merged output file does not exist: %s", output_fn)

  def gatk_preprocess_sample(self, sample, genome=None, libtype=None, wait=True):
    '''
    Convenience method to look up all the libraries for a sample and
    submit them to gatk_preprocess_libraries.
    '''
    libs = Library.objects.filter(sample__name=sample)
    if libtype is not None:
      libs = libs.filter(libtype__code__iexact=libtype)
    if libs.count() == 0:
      raise StandardError("No libraries found for sample %s (libtype %s)"
                          % (sample, libtype))

    libcodes = [ lib.code for lib in libs ]
    self.gatk_preprocess_libraries(libcodes, genome, wait=wait)

  def gatk_preprocess_libraries(self, libcodes, genome=None, wait=True):
    '''
    Runs our standard GATK preprocessing pipeline on a set of
    libraries. Sanity checks are made that the libraries all come from
    the same sample, are of the same type, and that the alignments are
    all against the same genome.
    '''
    libs = Library.objects.filter(code__in=libcodes)

    if libs.count() == 0:
      raise StandardError("No libraries match libcodes %s" % ",".join(libcodes))

    # Quick sanity check.
    indivs = list(set([ lib.sample.name for lib in libs]))
    if len(indivs) > 1:
      raise ValueError("Libraries come from multiple individual samples: %s"
                       % ", ".join(indivs))

    bams = Alnfile.objects.filter(alignment__lane__library__code__in=libcodes,
                                  filetype__code='bam')
    if genome is not None:
      bams = bams.filter(alignment__genome__code=genome)

    # Another sanity check.
    alngens = list(set([ bam.alignment.genome.code for bam in bams ]))
    if len(alngens) > 1:
      raise ValueError("Alignments found against multiple genomes: %s"
                       % ", ".join(alngens))

    # Another sanity check.
    libtypes = list(set([ bam.alignment.lane.library.libtype.code for bam in bams ]))
    if len(libtypes) > 1:
      raise ValueError("Alignments found against multiple library types: %s"
                       % ", ".join(libtypes))

    # Another sanity check.
    tissues = list(set([ bam.alignment.lane.library.sample.tissue.name for bam in bams ]))
    if len(tissues) > 1:
      raise ValueError("Alignments found against multiple tissues: %s"
                       % ", ".join(tissues))

    # And yet another sanity check.
    LOGGER.info("Validating bam file checksums.")
    for bam in bams:
      md5 = checksum_file(bam.repository_file_path, unzip=False)
      if md5 != bam.checksum:
        raise ValueError("Checksum for bam file %s does not agree with repository: (%s, %s)"
                         % (bam.filename, md5, bam.checksum))

    LOGGER.info("Count of %d bam files found for sample individual %s",
                bams.count(), indivs[0])
    merged_fn = "%s.bam" % (sanitize_samplename(indivs[0]),)
  
    # Now we merge the files.
    self.samtools_merge_bams([ bam.repository_file_path for bam in bams ],
                             merged_fn)

    self.gatk_preprocess_bam(merged_fn, bams[0].alignment, wait=wait)

  def gatk_preprocess_bam(self, merged_fn, aln=None, wait=True, genome=None):
    '''
    An alternative entry point which can be used if, say, one's
    pipeline has crashed and one doesn't want to wait for the bam file
    merge step to run again (given that this typically takes 7-8 hours
    for a full 40x WGS bam, it's sometimes worth doing this).
    '''
    # If not already supplied, retrieve one of the alignments; we
    # assume here that all alignments for the merged bam point to the
    # same library and alignment genome.
    if aln is None:
      with AlignmentFile(filename=merged_fn) as bamhandle:
        rgroups = bamhandle.header.get('RG', [])
      aln = retrieve_readgroup_alignment(rgroups[0], genome=genome)

    lib = aln.lane.library

    # Transfer the output to the cluster.
    cluster_merged = self.cluster_filename(merged_fn)
    LOGGER.info("Transfering files to cluster...")
    self.submitter.remote_copy_files(filenames=[ merged_fn ],
                                     destnames=[ cluster_merged ])

    finalpref = re.sub(' ', '_', ("%s%s_%s_"
                                  % (self.outprefix, lib.sample.tissue.name, lib.libtype.code)))
    (finalbam, finaljob) =\
        self.submit_cluster_jobs(cluster_merged,
                                 samplename=lib.sample.name,
                                 genobj=aln.genome,
                                 finalprefix=finalpref)

    # Check the expected output file is present in cwd.
    if wait:
      self.wait_on_cluster(jobs=[ finaljob ], finalbam=finalbam, localbam=merged_fn)

  def gatk_preprocess_free_bamfile(self, bamfile, samplename=None,
                                   genome=None, wait=True):
    '''
    Bare-bones processing of a bam file which is not associated with
    the repository. We assume as little as possible about this file.
    '''
    if samplename is None:
      samplename = os.path.splitext(bamfile)[0]

    genobj = None if genome is None else Genome.objects.get(code=genome)

    cluster_merged = self.cluster_filename(bamfile)
    LOGGER.info("Transfering files to cluster...")
    self.submitter.remote_copy_files(filenames=[ bamfile ],
                                     destnames=[ cluster_merged ])

    (finalbam, finaljob) = \
        self.submit_cluster_jobs(cluster_merged,
                                 samplename=samplename,
                                 genobj=genobj,
                                 finalprefix=self.outprefix)

    # Check the expected output file is present in cwd.
    if wait:
      self.wait_on_cluster(jobs=[ finaljob ], finalbam=finalbam, localbam=bamfile)

  def submit_markduplicates_job(self, cluster_merged):
    '''
    Submit a picard MarkDuplicates job to the cluster. This method
    will also submit a clean-up job to delete the MarkDuplicates log
    file.
    '''
    LOGGER.info("Submitting MarkDuplicates job")
    rootname = os.path.splitext(cluster_merged)[0]
    dupmark_fn  = "%s_dupmark.bam" % (rootname,)
    dupmark_log = "%s_dupmark.log" % (rootname,)
    dupmark_bai = "%s_dupmark.bai" % (rootname,)

    # The following limits java heap memory usage, which is obviously
    # important on the cluster. There are unconfirmed reports that
    # this actually helps with picard's stability as well.
    cmd = ('picard', '--Xmx', '8g', # requires picard python wrapper
           'MarkDuplicates',
           'I=%s' % quote(cluster_merged),
           'O=%s' % quote(dupmark_fn),
           'TMP_DIR=%s' % CONFIG.clusterworkdir,
           'VALIDATION_STRINGENCY=SILENT',
           'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512',
           'M=%s' % (quote(dupmark_log),))
    mdjob = self.submitter.submit_command(cmd=cmd,
                                          mem=10000,
                                          auto_requeue=False)

    # Cleanup job.
    cmd = ('rm', quote(dupmark_log))
    self.submitter.submit_command(cmd, depend_jobs=[ mdjob ])

    return (mdjob, dupmark_fn, dupmark_bai)

  def submit_buildbamindex_job(self, dupmark_fn, mdjob=None):
    '''
    Submit a picard BuildBamIndex job to the cluster.
    '''
    LOGGER.info("Submitting BuildBamIndex job")
    cmd = ('picard', '--Xmx', '8g',
           'BuildBamIndex',
           'VALIDATION_STRINGENCY=SILENT',
           'I=%s' % quote(dupmark_fn))
    depjobs = [] if mdjob is None else [ mdjob ]
    bijob = self.submitter.submit_command(cmd=cmd,
                                          mem=10000,
                                          depend_jobs=depjobs)

    return bijob

  def submit_gatk_pipeline_job(self, dupmark_fn, genobj, bijob=None, finalprefix=''):
    '''
    Submit a GATK pipeline job to the cluster. This method reads in
    the config.xml file in the top-level GATK pipeline installation
    directory, and writes a modified version to the cluster working
    area in which several variables are customised to this run. The
    cluster config file is deleted upon completion of the pipeline.

    In principle, this system can be used to run arbitrary pipeline
    components by altering the input config.xml; see the GATK pipeline
    documentation for details. Note that downstream code in this class
    may need to be altered depending on the outputs of any changed
    workflow (i.e., this class has expectations of the GATK output).
    '''
    LOGGER.info("Building GATK pipeline config file")
    tmpdir = os.path.join(CONFIG.clusterworkdir,
                          "%d_%d_gatk_tmp" % (os.getpid(), int(time())))
    conffile = self.create_instance_config(inputbam=dupmark_fn, # no quoting needed.
                                           tmpdir=tmpdir,
                                           outdir=CONFIG.gatk_cluster_output,
                                           reference=genobj.fasta_path,
                                           finalprefix=finalprefix)

    # Deliberately wait a couple of seconds so we definitively avoid
    # race conditions in tmpdir and conffile naming.
    sleep(2)

    # Then call gatk-pipeline/bin/run-pipeline --mode lsf
    # run_config.xml. Note that we need to set JAVA_HOME to the
    # correct location. Also, if we set the temp directory to a
    # pid-specific name then we can add --remove-temp to this to get a
    # better cleanup.
    LOGGER.info("Submitting GATK pipeline job")
    cmd = (os.path.join(CONFIG.gatk_cluster_root, 'bin', 'run-pipeline'),
           '--mode', 'lsf', '--remove-temp', conffile)
    environ = {'JAVA_HOME' : CONFIG.gatk_cluster_java_home}
    depjobs = [] if bijob is None else [ bijob ]
    gatkjob = self.submitter.submit_command(cmd, depend_jobs=depjobs,
                                            environ=environ)

    cmd = ('rm', conffile)
    self.submitter.submit_command(cmd, depend_jobs=[ gatkjob ])

    return gatkjob

  def submit_cluster_jobs(self, cluster_merged, samplename, genobj,
                          finalprefix):
    '''
    Submit the cluster jobs necessary to bridge between the samtools
    merge and the GATK pipeline invocation; also invokes the GATK
    pipeline, submits various output data transfer and cleanup jobs.

    Returns the name (but not the full path) of the local output bam file,
    and the number of the final cluster job which transferred the
    file back to the local host.
    '''
    # N.B. local var out_fns is assumed to always be in the order
    # (bam, bai).
    
    samplename = sanitize_samplename(samplename)

    if self.markduplicates:

      # Run MarkDuplicates
      (mdjob, dupmark_fn, dupmark_bai) =\
         self.submit_markduplicates_job(cluster_merged)

      # Cleanup job.
      cmd = ('rm', quote(cluster_merged))
      self.submitter.submit_command(cmd, depend_jobs=[ mdjob ])

      # Don't overwrite original input otherwise we have no intrinsic
      # control to test for MarkDuplicates failure.
      # Run BuildBamIndex
      lastjob     = self.submit_buildbamindex_job(dupmark_fn, mdjob)
      out_fns     = [ dupmark_fn, dupmark_bai ]
      gatkinput   = dupmark_fn

    else:

      # Just run BuildBamIndex so that the GATK pipeline is satisfied.
      lastjob     = self.submit_buildbamindex_job(cluster_merged)
      rootname    = os.path.splitext(cluster_merged)[0]
      cluster_bai = "%s.bai" % (rootname,)
      out_fns     = [ cluster_merged, cluster_bai ]
      gatkinput   = cluster_merged

    if self.gatkpipe:

      if genobj is None:
        raise ValueError("Genome object not specified in GATK pipeline call.")

      # Run the GATK pipeline.
      lastjob = self.submit_gatk_pipeline_job(gatkinput,
                                              genobj=genobj,
                                              bijob=lastjob,
                                              finalprefix=finalprefix)
      # Cleanup job.
      cmd = 'rm %s' % " ".join([quote(x) for x in out_fns])
      self.submitter.submit_command(cmd, depend_jobs=[ lastjob ])

      out_fns = [ os.path.join(CONFIG.gatk_cluster_output,
                               "%s%s.%s" % (finalprefix, samplename, suffix,))
                  for suffix in ('bam','bai') ]

    # Copy the output back to local cwd. Also cleanup, but only if
    # transfer was successful.
    LOGGER.info("Submitting output bam file transfer job")
    cmd = " && ".join([ \
      self.return_file_to_localhost(fname,
                                    _make_local_bamfile_name(fname, samplename, finalprefix),
                                    donefile=True,
                                    execute=False) for fname in out_fns ])
    cmd += ' && rm %s' % " ".join([quote(x) for x in out_fns])
    sshjob = self.submitter.submit_command(cmd, depend_jobs=[ lastjob ])

    return (_make_local_bamfile_name(out_fns[0], samplename, finalprefix), sshjob)

  def create_instance_config(self, inputbam, tmpdir, outdir,
                             reference, finalprefix=''):
    '''
    Read in the edited template config file from the GATK pipeline,
    modify a couple of run-specific variables, and write it out to the
    GATK cluster input directory.

    Run-specific variables:
    1. genome fasta path
    2. input bamFiles
    3. temp directory name
    4. (to be completed TODO) various known indel/snv vcfs
    '''

    # TODO consider using the original config template as provided by
    # the pipeline package, rather than our custom-edited version here.
    conffile = os.path.join(CONFIG.gatk_cluster_root, 'config.xml')
    new_conffile = os.path.join(CONFIG.gatk_cluster_input,
                                "%d_%d_config.xml" % (os.getpid(), int(time())))

    def_conf = ET.parse(conffile)

    var_elem = def_conf.find('.//variables')

    ref_elem = var_elem.find('./referenceSequence')
    ref_elem.text = str(reference)
    
    dir_elem = var_elem.find('./outputPrefix')
    dir_elem.text = str(finalprefix)

    dir_elem = var_elem.find('./outputDir')
    dir_elem.text = str(outdir)

    bam_elem = var_elem.find('./bamFiles')
    bam_elem.text = str(inputbam)

    dir_elem = var_elem.find('./work')
    dir_elem.text = str(tmpdir)

    tmpconf = NamedTemporaryFile(delete=False, dir=CONFIG.tmpdir)
    def_conf.write(tmpconf.name, encoding='ISO-8859-1', xml_declaration=True)
    self.submitter.remote_copy_files(filenames=[ tmpconf.name ],
                                     destnames=[ new_conffile ])
    os.unlink(tmpconf.name)

    return new_conffile

