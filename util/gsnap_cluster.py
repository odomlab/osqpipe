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
A script which can be used to run GSNAP on the cluster over a set of
paired sequencing fastq files.
'''

import re
import os
import gzip
from osqpipe.pipeline.bwa_runner import ClusterJobManager
from osqutil.utilities import is_zipped
from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

################################################################################

class GsnapManager(ClusterJobManager):
  '''
  Class used to manage the submission of gsnap jobs to the cluster.
  '''
  def __init__(self, outdir='gsnap', *args, **kwargs):

    super(GsnapManager, self).__init__(*args, **kwargs)

    self.outdir = outdir

  @staticmethod
  def count_fastq_reads(fqfile, file_zipped=False):
    '''
    Count the number of reads in a fastq file.
    '''
    if file_zipped:
      fqfh = gzip.open(fqfile, 'rb')
    else:
      fqfh = open(fqfile, 'r')

    count = 0
    for _line in fqfh:
      count +=1

    fqfh.close()
    return int(count / 4)

  def submit_gsnap(self, files, genome, indexdir, number, queue, wait=False):
    '''
    Submit gsnap alignment jobs, samtools merge and cleanup jobs to the cluster.
    '''
    assert len(files) in (1,2)

    keeptype = 'concordant' if len(files) == 2 else 'unpaired'

    bjobs = []

    # Only test the zippedness of our inputs once.
    use_gunzip = True if is_zipped(files[0]) else False

    if number is None:
      LOGGER.info("Counting reads in fastq input files...")
      number = int((self.count_fastq_reads(files[0], use_gunzip) / 1e6) * 20)
      LOGGER.info("Will start %d jobs on cluster.", number)

    libmatch = re.match(r'^(do\d+)_', os.path.basename(files[0]))
    if libmatch:
      libcode = libmatch.group(1)
    else:
      raise ValueError("Unable to parse library code from file name %s"
                       % files[0])

    # Copy files to cluster workdir.
    clfiles = [ '%s_%s' % (self.namespace, os.path.basename(x)) for x in files ]
    for filenum in range(len(files)):
      if not self.cluster_file_exists(clfiles[filenum]):
        LOGGER.info("Copying fastq file to cluster working directory...")
        self.submitter.remote_copy_files([files[filenum]], [clfiles[filenum]])

    # Launch the gsnap job array.
    prefix = (r'%s/%s_%s.\\\$((\\\$LSB_JOBINDEX - 1))'
              % (self.outdir, self.namespace, libcode))

    # Create the working directory, if it's not already there.
    cmd  = r'mkdir -p %s' % self.outdir

    # Build the gsnap command. Note that the '-B' option can be
    # switched to 5 if there's too much disk IO; this has the
    # downside of increasing memory requirements though.
    cmd += ((r' && gsnap --part=\\\$((\\\$LSB_JOBINDEX - 1))/%d -d %s -D %s'
             % (number, genome, indexdir))
            + r' -m 2 --trim-mismatch-score=0 -N 1 -n 1'
            + r' -E 100 --pairexpect=160 -B 4'
            + r' --antistranded-penalty=1 --quality-protocol=sanger'
            + r' --nofails -A sam'
            + (r' --split-output=%s --nthreads=8' % prefix))

    if use_gunzip:
      cmd += ' --gunzip'

    cmd += " " + " ".join(clfiles)

    # The samtools view (sam->bam) step is actually where we assume
    # paired-end sequencing.
    for tag in ('uniq', 'mult'):
      samfile = r'%s.%s_%s' % (prefix, keeptype, tag)

      # Sort here rather than later as it's more efficient.
      cmd += (r' && samtools view -b -S -o - %s | samtools sort -o %s.bam -'
              % (samfile, samfile))

    # Clean up the sam files.
    cmd += r' && rm %s.nomapping' % prefix
    cmd += ((r' && rm %s.' % prefix)
            + r'{unpaired,paired,halfmapping,concordant}_'
            + r'{transloc,mult,circular,uniq,uniq_scr,uniq_long,uniq_inv,uniq_circular}')

    cmd = (r'if [ ! -f %s.%s_mult.bam ]; then %s; else echo Skipping %s; fi'
           % (prefix, keeptype, cmd, prefix))

    LOGGER.info("Submitting job array of %d jobs.", number)
    jobname = 'gsnap[1-%d]%%%d' % (number, self.throttle)
    LOGGER.debug("Job name: %s", jobname)
    LOGGER.debug("Command: %s", cmd)
    bjobs.append(self.submitter.submit_command(cmd, mem=self.memsize * 1024,
                                               queue=queue,
                                               auto_requeue=False,
                                               jobname=jobname))

    mergejob = self._merge_bams(bjobs, libcode, queue,
                                njobs=number, keeptype=keeptype)

    fastq_cleanup = r'rm %s' % " ".join(clfiles)

    if wait is True:
      self.wait_on_cluster([mergejob], fastq_cleanup)
      return

    else:
      LOGGER.info("Submitting fastq cleanup job...")
      jobid = self.submitter.submit_command(fastq_cleanup,
                                            queue=queue,
                                            depend_jobs=[mergejob])
      return jobid

  def _merge_bams(self, bjobs, libcode, queue,
                  limit=500, njobs=None, keeptype='concordant'):
    '''
    Merge the output bams from a set of cluster jobs.
    '''
    # If we're restarting a partial run, len(bjobs) might be smaller
    # than we'd like.
    if njobs is None:
      LOGGER.warn("Bamfile merge step will take number of jobs from"
                  + " the LSF jobid list.")
      njobs  = len(bjobs)

    prefix = '%s/%s_%s' % (self.outdir, self.namespace, libcode)
    suffix = '%s_{uniq,mult}.bam' % keeptype

    # Check that the number of bam files found equals njobs*2.
    bampatt = '%s.{0..%d}.%s' % (prefix, njobs - 1, suffix)
    precmd  = (r'if [ \\\`ls %s | wc -l\\\` == %d ]; then '
               % (bampatt, njobs * 2))
    postcmd = r'; else echo Unsufficient bam files created. Stopping.; fi'

    # If we have to merge more than ~1,000 files, the OS will complain
    # due to default ulimit settings. We merge such file sets in
    # stages. I'm assuming here that we won't ever need to merge
    # 1,000,000 files; otherwise, this method will need some kind of recursion.
    div = 1
    while njobs / div > limit:
      div += 1

    mergejobs = []
    batchsize = njobs / div
    for batchnum in range(div):

      start = batchnum * batchsize

      # Safety check to make sure our integer math doesn't drop stray
      # runs in the last batch.
      end   = njobs-1 if batchnum == div-1 else ((batchnum+1) * batchsize) - 1

      filepatt = '%s.{%d..%d}.%s' % (prefix, start, end, suffix)

      # Intermediary bam files in the form "$prefix.0.bam, $prefix.1.bam,..."
      if div == 1:
        output = '%s.bam' % (prefix,)
      else:
        output = '%s.%d.bam' % (prefix, batchnum)

      # We assume here that the bam files are already sorted.
      cmd = ('samtools merge %s %s' % (output, filepatt))
      cmd = precmd + cmd + postcmd

      LOGGER.info("Submitting bamfile merge job (%d/%d)...", batchnum+1, div)
      LOGGER.debug(cmd)
      jobid = self.submitter.submit_command(cmd, mem=20000,
                                            queue=queue,
                                            auto_requeue=False,
                                            depend_jobs=bjobs)
      mergejobs.append(jobid)

    # Submit a job to merge the intermediary bams into a final output bam file.
    if len(mergejobs) > 1:
      filepatt = '%s.{0..%d}.bam' % (prefix, len(mergejobs) - 1)

      # Again, we assume the bam files are sorted.
      cmd = ('samtools merge %s.bam %s && rm %s'
             % (prefix, filepatt, filepatt))
      cmd = precmd + cmd + postcmd

      LOGGER.info("Submitting final bamfile merge job...")
      LOGGER.debug(cmd)
      jobid = self.submitter.submit_command(cmd, mem=20000,
                                            queue=queue,
                                            auto_requeue=False,
                                            depend_jobs=mergejobs)
    else:
      jobid = mergejobs[0]

    # Final cleanup of initial bam files.
    cmd = precmd + 'rm %s' % bampatt + postcmd
    LOGGER.info("Submitting bamfile cleanup job...")
    LOGGER.debug(cmd)
    jobid = self.submitter.submit_command(cmd,
                                          queue=queue,
                                          auto_requeue=False,
                                          depend_jobs = [jobid])

    return jobid

################################################################################

if __name__ == '__main__':

  from argparse import ArgumentParser

  P = ArgumentParser(description='Wrapper script to copy fastq files to the'
                     + ' cluster, run GSNAP over them, and to combine'
                     + ' the outputs.')

  P.add_argument('files', metavar='<filenames>', type=str, nargs='+',
                      help='The name of the fastq file(s) to process.\n These'
                      + ' must be either one (SE) or two (PE) fastq files.')

  P.add_argument('-d', '--genome', dest='genome', type=str, required=True,
                 help='The genome code (as understood by GMAP/GSNAP) to align'
                 + ' against.')

  P.add_argument('-D', '--indexdir', dest='indexdir', type=str, required=True,
                 help='The cluster directory containing GSNAP genome indices.')

  P.add_argument('-n', '--number', dest='number', type=int,
                 help='The number of cluster jobs to run. If omitted, a'
                 + ' best-guess number will be derived from the input'
                 + ' fastq files.')

  P.add_argument('--namespace', dest='namespace', type=str,
                 default=str(os.getpid()),
                 help='A namespace argument, used as a prefix for filenames'
                 + ' on the cluster. Defaults to using the script process'
                 + ' ID (pid).')

  P.add_argument('--throttle', dest='throttle', type=int, default=40,
                 help='The number of concurrent gsnap jobs which may be running'
                 + ' on the cluster (Default: 40).')

  P.add_argument('--memory', dest='memsize', type=int, default=20,
                 help='The memory (in GB) to request for each gsnap job.'
                 + ' Default is 20GB.')

  P.add_argument('--queue', dest='queue', type=str,
                 help='The name of the cluster queue to submit to.')

  P.add_argument('-i', '--ssh-key', dest='ssh_key', type=str,
                 help='The path to a (passwordless) ssh key used by the cluster'
                 + ' account to communicate job completion to our local host.')

  P.add_argument('--wait', dest='wait', action='store_true',
                 help='Whether to wait for the final cluster job to complete'
                 + ' before exiting the script. Useful in shell loops.')

  ARGS = P.parse_args()

  GMAN = GsnapManager(namespace = ARGS.namespace,
                      throttle  = ARGS.throttle,
                      memsize   = ARGS.memsize,
                      ssh_key   = ARGS.ssh_key)

  GMAN.submit_gsnap(files    = ARGS.files,
                    genome   = ARGS.genome,
                    indexdir = ARGS.indexdir,
                    number   = ARGS.number,
                    queue    = ARGS.queue,
                    wait     = ARGS.wait)
