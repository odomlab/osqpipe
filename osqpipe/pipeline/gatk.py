'''
Code used to manage the GATK cleanup parts of our HCC pipeline.
'''

import os
import re
from subprocess import Popen, PIPE, STDOUT, CalledProcessError
from django.db import transaction
from osqpipe.models import Alnfile, Library, Alignment
from osqpipe.pipeline.samtools import count_bam_reads
from osqpipe.pipeline.utilities import call_subprocess, checksum_file, \
    sanitize_samplename
from pipes import quote
from osqpipe.pipeline.config import Config
from osqpipe.pipeline.bwa_runner import ClusterJobManager

import xml.etree.ElementTree as ET
from tempfile import NamedTemporaryFile
from pysam import AlignmentFile
from shutil import move

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging('gatk')
CONFIG = Config()

OUTPREF='IR_BQSR_'

################################################################################
# A handful of utility functions.
def retrieve_readgroup_alignment(rgroup, genome=None):
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

  if alns.count() > 1:
    raise ValueError("Multiple Alignments match read group and genome parameters.")
  elif alns.count() == 0:
    raise StandardError("No Alignments found to match read group and genome parameters.")
  else:
    return alns[0]

def check_bam_readcount(bam, maln):
  '''
  In principle, the total reads returned by the GATK pipeline should
  be the sum of the reads in the original fastq files. We check that
  here. Arguments are: bam (the GATK output bam file), maln (the
  MergedAlignment db object).
  '''
  expected = sum([ aln.lane.total_passedpf for aln in maln.alignments.all() ])
  numreads = count_bam_reads(bam)

  ## See how things pan out: if we have to relax this check, here
  ## would be a good place to start (i.e., raise a warning rather than
  ## an Exception).
  if numreads != expected:
    message = ("Number of reads in bam file is differs from that in "
               + "fastq file: %d (bam) vs %d (fastq)")
    raise ValueError(message % (numreads, expected))

@transaction.commit_on_success
def update_bam_readgroups(bam):
  '''
  A lightweight method for updating bam read group information to
  match the current repository annotation.
  '''
  # Support both the filename or an Alnfile object being
  # passed. MergedAlnfile isn't really appropriate here.
  if issubclass(type(bam), Alnfile):
    bam = Alnfile.objects.get(id=bam.id)
  else:
    bam = Alnfile.objects.get(filename=bam)

  if bam.filetype.code != 'bam':
    raise ValueError("Function requires bam file, not %s (%s)"
                     % (bam.filetype.code, bam.filename))

  library = bam.alignment.lane.library

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
      fields['PU'] = bam.alignment.lane.lanenum
      fields['LB'] = library.code
      fields['SM'] = sanitize_samplename(library.individual)
      fields['CN'] = bam.alignment.lane.facility.code
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
  def cluster_filename(self, fname):
    '''
    Given a local filename, return a file path suitable for use on the
    cluster.
    '''
    fname  = "%d_%s" % (os.getpid(), fname)
    clpath = os.path.join(CONFIG.gatk_cluster_input, fname)
    return clpath

  def samtools_merge_bams(self, bams, output_fn):
    '''
    Use samtools to merge a set of bams locally.
    '''
    LOGGER.info("Using samtools to merge bam files: %s",
                ", ".join([ os.path.basename(bam) for bam in bams ]))

    if len(bams) == 0:
      raise ValueError("Zero input bam files for merging.")
  
    # We assume our input bam files are appropriately sorted (which,
    # coming from the repository, they should be).
    cmd = ['samtools', 'merge', output_fn] + bams

    call_subprocess(cmd, path=CONFIG.hostpath)

    if not os.path.exists(output_fn):
      raise StandardError("Merged output file does not exist: %s", output_fn)

  def gatk_preprocess_libraries(self, libcodes, genome=None,
                                outprefix=OUTPREF):
    '''
    Main entry point for the class.
    '''
    libs = Library.objects.filter(code__in=libcodes)

    # Quick sanity check.
    indivs = list(set([ lib.individual for lib in libs]))
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
    tissues = list(set([ bam.alignment.lane.library.tissue.name for bam in bams ]))
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

    self.gatk_preprocess_bam(merged_fn, outprefix, bams[0].alignment)

  def gatk_preprocess_bam(self, merged_fn, outprefix=OUTPREF, aln=None, wait=True):
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
      aln = retrieve_readgroup_alignment(rgroups[0])

    lib = aln.lane.library

    # Transfer the output to the cluster.
    cluster_merged = self.cluster_filename(merged_fn)
    LOGGER.info("Transfering files to cluster...")
    self.submitter.remote_copy_files(filenames=[ merged_fn ],
                                     destnames=[ cluster_merged ])

    finalpref = re.sub(' ', '_', ("%s%s_%s_"
                                  % (outprefix, lib.tissue.name, lib.libtype.code)))
    (finalbam, finaljob) =\
        self.submit_cluster_jobs(cluster_merged,
                                 samplename=lib.individual,
                                 genobj=aln.genome,
                                 outprefix=finalpref)

    # Check the expected output file is present in cwd.
    if wait:
      self.wait_on_cluster([ finaljob ])
      finaldone = "%s.done" % finalbam
      if not (os.path.exists(finalbam) and os.path.exists(finaldone)):
        raise StandardError("Expected output file %s not found, or %s not created."
                            % (finalbam, finaldone))

      # Delete local files (only if we're waiting on cluster
      # though). Note that if the cluster jobs fail, the local merged
      # bam file should not be deleted.
      LOGGER.info("Deleting local merged bam file")
      os.unlink(merged_fn)

  def submit_markduplicates_job(self, cluster_merged):
    '''
    Submit a picard MarkDuplicates job to the cluster. This method
    will also submit a clean-up job to delete the MarkDuplicates log
    file.
    '''
    LOGGER.info("Submitting MarkDuplicates job")
    cluster_merged = cluster_merged
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

  def submit_buildbamindex_job(self, dupmark_fn, mdjob):
    '''
    Submit a picard BuildBamIndex job to the cluster.
    '''
    LOGGER.info("Submitting BuildBamIndex job")
    cmd = ('picard', '--Xmx', '8g',
           'BuildBamIndex',
           'VALIDATION_STRINGENCY=SILENT',
           'I=%s' % quote(dupmark_fn))
    bijob = self.submitter.submit_command(cmd=cmd,
                                          mem=10000,
                                          depend_jobs=[ mdjob ])

    return bijob

  def submit_gatk_pipeline_job(self, dupmark_fn, bijob, genobj,
                               outprefix=OUTPREF):
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
    tmpdir = os.path.join(CONFIG.clusterworkdir, "%d_gatk_tmp" % os.getpid())
    conffile = self.create_instance_config(inputbam=dupmark_fn, # no quoting needed.
                                           tmpdir=tmpdir,
                                           outdir=CONFIG.gatk_cluster_output,
                                           reference=genobj.fasta_path,
                                           outprefix=outprefix)

    # Then call gatk-pipeline/bin/run-pipeline --mode lsf
    # run_config.xml. Note that we need to set JAVA_HOME to the
    # correct location. Also, if we set the temp directory to a
    # pid-specific name then we can add --remove-temp to this to get a
    # better cleanup.
    LOGGER.info("Submitting GATK pipeline job")
    cmd = (os.path.join(CONFIG.gatk_cluster_root, 'bin', 'run-pipeline'),
           '--mode', 'lsf', '--remove-temp', conffile)
    environ = {'JAVA_HOME' : CONFIG.gatk_cluster_java_home}
    gatkjob = self.submitter.submit_command(cmd, depend_jobs=[ bijob ],
                                            environ=environ)

    cmd = ('rm', conffile)
    self.submitter.submit_command(cmd, depend_jobs=[ gatkjob ])

    return gatkjob

  def submit_cluster_jobs(self, cluster_merged, samplename, genobj,
                          outprefix=OUTPREF):
    '''
    Submit the cluster jobs necessary to bridge between the samtools
    merge and the GATK pipeline invocation; also invokes the GATK
    pipeline, submits various output data transfer and cleanup jobs.
    '''
    # Run MarkDuplicates
    (mdjob, dupmark_fn, dupmark_bai) =\
        self.submit_markduplicates_job(cluster_merged)

    # Cleanup job.
    cmd = ('rm', quote(cluster_merged))
    self.submitter.submit_command(cmd, depend_jobs=[ mdjob ])

    # Don't overwrite original input otherwise we have no intrinsic
    # control to test for MarkDuplicates failure.
    # Run BuildBamIndex
    bijob = self.submit_buildbamindex_job(dupmark_fn, mdjob)

    # Run the GATK pipeline.
    gatkjob = self.submit_gatk_pipeline_job(dupmark_fn, bijob,
                                            genobj=genobj,
                                            outprefix=outprefix)

    # Cleanup job.
    cmd = ('rm', quote(dupmark_fn), quote(dupmark_bai))
    self.submitter.submit_command(cmd, depend_jobs=[ gatkjob ])

    # Copy the output back to local cwd.
    LOGGER.info("Submitting output bam file transfer job")
    finalbam = "%s%s.bam" % (outprefix, samplename,)
    clusterout = os.path.join(CONFIG.gatk_cluster_output,
                              "%s%s.bam" % (outprefix, samplename,))
    clusterbai = os.path.join(CONFIG.gatk_cluster_output,
                              "%s%s.bai" % (outprefix, samplename,))
    cmd = self.return_file_to_localhost(clusterout,
                                        finalbam,
                                        donefile=True,
                                        execute=False)
    sshjob = self.submitter.submit_command(cmd, depend_jobs=[ gatkjob ])

    # Cleanup job.
    cmd = ('rm', quote(clusterout), quote(clusterbai))
    self.submitter.submit_command(cmd, depend_jobs=[ sshjob ])

    return (finalbam, sshjob)

  def create_instance_config(self, inputbam, tmpdir, outdir,
                             reference, outprefix=OUTPREF):
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
                                "%d_config.xml" % os.getpid())

    def_conf = ET.parse(conffile)

    var_elem = def_conf.find('.//variables')

    ref_elem = var_elem.find('./referenceSequence')
    ref_elem.text = str(reference)
    
    dir_elem = var_elem.find('./outputPrefix')
    dir_elem.text = str(outprefix)

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

