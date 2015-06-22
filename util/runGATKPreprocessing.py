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

  def cluster_filename(self, fname):
    fname  = "%d_%s" % (os.getpid(), fname)
    clpath = os.path.join(CONFIG.gatk_cluster_input, fname)
    return clpath

  def samtools_merge_bams(self, bams, output_fn):

    LOGGER.info("Using samtools to merge bam files: %s", ", ".join([ os.path.basename(bam) for bam in bams ]))
  
    # We assume our input bam files are appropriately sorted (which,
    # coming from the repository, they should be).
    cmd = ['samtools', 'merge', output_fn] + bams

    call_subprocess(cmd, path=CONFIG.hostpath)

  def gatk_preprocess_libraries(self, libcodes, genome=None, outprefix='IR_BQSR_'):
    '''
    Main entry point.
    '''
    libs   = Library.objects.filter(code__in=libcodes)

    # Quick sanity check.
    indivs = list(set([ lib.individual for lib in libs])) # TODO count() somehow?
    if len(indivs) > 1:
      raise ValueError("Libraries come from multiple individual samples: %s" % ", ".join(indivs))

    bams = Alnfile.objects.filter(alignment__lane__library__code__in=libcodes, filetype__code='bam')
    if genome is not None:
      bams = bams.filter(alignment__genome__code=genome)

    # Another sanity check.
    alngens = list(set([ bam.alignment.genome.code for bam in bams ])) # TODO count() somehow?
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
    output_fn = "%s.bam" % (indivs[0],)
  
    # Now we merge the files.
#    self.samtools_merge_bams([ bam.repository_file_path for bam in bams ], output_fn)

    if not os.path.exists(output_fn):
      raise StandardError("Merged output file does not exist: %s", output_fn)

    # Transfer the output to the cluster.
#    cluster_input = self.cluster_filename(output_fn)
    cluster_input = '/lustre/dolab/fnc-odompipe/pipeline/gatk_input/5395_80400.bam'
    LOGGER.info("Transfering files to cluster...")
#    self.submitter.remote_copy_files(filenames=[ output_fn ],
#                                     destnames=[ cluster_input ])

    # Run MarkDuplicates
    LOGGER.info("Submitting MarkDuplicates job")
    rootname = os.path.splitext(cluster_input)[0]
    md_output = "%s_dupmark.bam" % (rootname,)
    md_log    = "%s_dupmark.log" % (rootname,)
    cmd = ('picard',
           'MarkDuplicates',
           'I=%s' % cluster_input,
           'O=%s' % md_output,
           'TMP_DIR=%s' % CONFIG.clusterworkdir,
           'VALIDATION_STRINGENCY=SILENT',
           'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512',
           'M=%s' % (md_log,))
#    mdjob = self.submitter.submit_command(cmd=cmd,
#                                          mem=10000,
#                                          auto_requeue=False)

    # Don't overwrite original input otherwise we have no intrinsic
    # control to test for MarkDuplicates failure.
    # Run BuildBamIndex
    cmd = ('picard',
           'BuildBamIndex',
           'I=%s' % md_output)
#    bijob = self.submitter.submit_command(cmd=cmd,
#                                          depend_jobs=[ mdjob ])

    # Run the GATK pipeline.
    genobj   = bams[0].alignment.genome
    conffile = self.create_instance_config(inputbam=md_output,
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
    cmd = (os.path.join(CONFIG.gatk_cluster_root, 'bin', 'run-pipeline'),
           '--mode', 'lsf', '--remove-temp', conffile)
    gatkjob = self.submitter.submit_command(cmd, # FIXME delete comment to reinstate dependency: depend_jobs=[ bijob ],
                                            environ={'JAVA_HOME' : CONFIG.gatk_cluster_java_home})

    # Copy the output back to cwd.
    finalbam = "%s.bam" % (indivs[0],)
    cmd = self.return_file_to_localhost("%s%s.bam" % (outprefix, indivs[0],),
                                        finalbam,
                                        donefile=True,
                                        execute=False)
    sshjob = self.submitter.submit_command(cmd, depend_jobs=[ gatkjob ])

    # Check the expected output file is present in cwd.
    self.wait_on_cluster([ sshjob ])
    finaldone = "%s.done" % finalbam
    if not os.path.exists(finalbam) and os.path.exists(finaldone):
      raise StandardError("Expected output file %s not found, or %s not created."
                          % (finalbam, finaldone))

    # Don't delete anything unless we're waiting on the cluster.
#    cmd = ('rm', cluster_input)
#    self.submitter.submit_command(cmd, depend_jobs=[ mdjob ])
#    cmd = ('rm', md_output, md_log)
#    self.submitter.submit_command(cmd, depend_jobs=[ gatkjob ])
#    cmd = ('rm', output_fn, finaldone, md_log, cluster_input)
#    self.submitter.submit_command(cmd, depend_jobs=[ sshjob ])

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
    # the pipeline package.
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
