#!/usr/bin/env python

'''
Script which bridges the gap between our standard sequencing pipeline
and the Bioinformatics Core GATK preprocessing pipeline.

The steps included in this script are:

1. Merge all the files from the specified libraries (which should all
come from the same HCC nodule).

2. Transfer the merged bam file to the cluster.

2. Run picard MarkDuplicates on the merged bam file.

3. Run picard BuildBamIndex on the output.

4. Start the GATK IndelRealigner-BaseRecalibrator pipeline as provided
by the Bioinformatics Core pipeline.

5. Submit cleanup jobs to transfer the output back to local host, and
delete working files on the cluster. The GATK log files are currently
retained in the working directory.
'''

import os
from osqpipe.models import Alnfile, Library
from osqpipe.pipeline.utilities import call_subprocess, checksum_file
from osqpipe.pipeline.config import Config
from osqpipe.pipeline.bwa_runner import ClusterJobManager

import xml.etree.ElementTree as ET
from tempfile import NamedTemporaryFile

from logging import INFO
from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)
CONFIG = Config()

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
  
    # We assume our input bam files are appropriately sorted (which,
    # coming from the repository, they should be).
    cmd = ['samtools', 'merge', output_fn] + bams

    call_subprocess(cmd, path=CONFIG.hostpath)

    if not os.path.exists(output_fn):
      raise StandardError("Merged output file does not exist: %s", output_fn)

  def gatk_preprocess_libraries(self, libcodes, genome=None, outprefix='IR_BQSR_'):
    '''
    Main entry point for the class.
    '''
    libs = Library.objects.filter(code__in=libcodes)

    # Quick sanity check.
    indivs = list(set([ lib.individual for lib in libs]))
    if len(indivs) > 1:
      raise ValueError("Libraries come from multiple individual samples: %s" % ", ".join(indivs))

    bams = Alnfile.objects.filter(alignment__lane__library__code__in=libcodes, filetype__code='bam')
    if genome is not None:
      bams = bams.filter(alignment__genome__code=genome)

    # Another sanity check.
    alngens = list(set([ bam.alignment.genome.code for bam in bams ]))
    if len(alngens) > 1:
      raise ValueError("Alignments found against multiple genomes: %s" % ", ".join(alngens))

    # And yet another sanity check.
    LOGGER.info("Validating bam file checksums.")
    for bam in bams:
      md5 = checksum_file(bam.repository_file_path, unzip=False)
      if md5 != bam.checksum:
        raise ValueError("Checksum for bam file %s does not agree with repository: (%s, %s)"
                         % (bam.filename, md5, bam.checksum))

    LOGGER.info("Count of %d bam files found for sample individual %s", bams.count(), indivs[0])
    merged_fn = "%s.bam" % (indivs[0],)
  
    # Now we merge the files.
    self.samtools_merge_bams([ bam.repository_file_path for bam in bams ], merged_fn)

    # Transfer the output to the cluster.
    cluster_merged = self.cluster_filename(merged_fn)
    LOGGER.info("Transfering files to cluster...")
    self.submitter.remote_copy_files(filenames=[ merged_fn ],
                                     destnames=[ cluster_merged ])

    genobj = bams[0].alignment.genome
    (finalbam, finaljob) = self.submit_cluster_jobs(cluster_merged,
                                                    samplename=indivs[0],
                                                    genobj=genobj,
                                                    outprefix=outprefix)

    # Check the expected output file is present in cwd.
    self.wait_on_cluster([ finaljob ])
    finaldone = "%s.done" % finalbam
    if not (os.path.exists(finalbam) and os.path.exists(finaldone)):
      raise StandardError("Expected output file %s not found, or %s not created."
                          % (finalbam, finaldone))

    # Delete local files (only if we're waiting on cluster though).
    LOGGER.info("Deleting local merged bam file")
    os.unlink(merged_fn)

  def submit_markduplicates_job(self, cluster_merged):

    LOGGER.info("Submitting MarkDuplicates job")
    rootname = os.path.splitext(cluster_merged)[0]
    dupmark_fn  = "%s_dupmark.bam" % (rootname,)
    dupmark_log = "%s_dupmark.log" % (rootname,)
    dupmark_bai = "%s_dupmark.bai" % (rootname,)
    cmd = ('picard',
           'MarkDuplicates',
           'I=%s' % cluster_merged,
           'O=%s' % dupmark_fn,
           'TMP_DIR=%s' % CONFIG.clusterworkdir,
           'VALIDATION_STRINGENCY=SILENT',
           'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512',
           'M=%s' % (dupmark_log,))
    mdjob = self.submitter.submit_command(cmd=cmd,
                                          mem=10000,
                                          auto_requeue=False)

    # Cleanup job.
    cmd = ('rm', dupmark_log)
    self.submitter.submit_command(cmd, depend_jobs=[ mdjob ])

    return (mdjob, dupmark_fn, dupmark_bai)

  def submit_buildbamindex_job(self, dupmark_fn, mdjob):

    LOGGER.info("Submitting BuildBamIndex job")
    cmd = ('picard',
           'BuildBamIndex',
           'I=%s' % dupmark_fn)
    bijob = self.submitter.submit_command(cmd=cmd,
                                          depend_jobs=[ mdjob ])

    return bijob

  def submit_gatk_pipeline_job(self, dupmark_fn, bijob, genobj, outprefix='IR_BQSR_'):

    LOGGER.info("Building GATK pipeline config file")
    conffile = self.create_instance_config(inputbam=dupmark_fn,
                                           tmpdir=os.path.join(CONFIG.clusterworkdir,
                                                               "%d_gatk_tmp" % os.getpid()),
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
    gatkjob = self.submitter.submit_command(cmd, depend_jobs=[ bijob ],
                                            environ={'JAVA_HOME' : CONFIG.gatk_cluster_java_home})

    cmd = ('rm', conffile)
    self.submitter.submit_command(cmd, depend_jobs=[ gatkjob ])

    return gatkjob

  def submit_cluster_jobs(self, cluster_merged, samplename, genobj, outprefix='IR_BQSR_'):

    # Run MarkDuplicates
    (mdjob, dupmark_fn, dupmark_bai) = self.submit_markduplicates_job(cluster_merged)

    # Cleanup job.
    cmd = ('rm', cluster_merged)
    self.submitter.submit_command(cmd, depend_jobs=[ mdjob ])

    # Don't overwrite original input otherwise we have no intrinsic
    # control to test for MarkDuplicates failure.
    # Run BuildBamIndex
    bijob = self.submit_buildbamindex_job(dupmark_fn, mdjob)

    # Run the GATK pipeline.
    gatkjob = self.submit_gatk_pipeline_job(dupmark_fn, bijob, genobj=genobj, outprefix=outprefix)

    # Cleanup job.
    cmd = ('rm', dupmark_fn, dupmark_bai)
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
    cmd = ('rm', clusterout, clusterbai)
    self.submitter.submit_command(cmd, depend_jobs=[ sshjob ])

    return (finalbam, sshjob)

  def create_instance_config(self, inputbam, tmpdir, outdir, reference, outprefix='IR_BQSR_'):
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
    ref_elem.text = str(reference) # Or unicode? encoding is not utf though
    
    dir_elem = var_elem.find('./outputPrefix')
    dir_elem.text = outprefix

    dir_elem = var_elem.find('./outputDir')
    dir_elem.text = outdir

    bam_elem = var_elem.find('./bamFiles')
    bam_elem.text = inputbam

    dir_elem = var_elem.find('./work')
    dir_elem.text = tmpdir

    tmpconf = NamedTemporaryFile(delete=False, dir=CONFIG.tmpdir)
    def_conf.write(tmpconf.name, encoding='ISO-8859-1', xml_declaration=True)
    self.submitter.remote_copy_files(filenames=[ tmpconf.name ],
                                     destnames=[ new_conffile ])
    os.unlink(tmpconf.name)

    return new_conffile

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description=\
               'Script to initiate the HCC GATK preprocessing pipeline.')
  
  PARSER.add_argument('libraries', metavar='<libcodes>', type=str, nargs='*',
                      help='The names of the libraries to merge and load into the pipeline.'
                      + ' All the files on the command line should come from the same HCC nodule.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=False,
                      help='The alignment genome used to filter the input files.')

  ARGS = PARSER.parse_args()

  PROC = GATKPreprocessor()
  PROC.gatk_preprocess_libraries(ARGS.libraries, genome=ARGS.genome)
