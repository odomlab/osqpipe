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
Script to handle the transfer of tumour and normal bam and bai
files to the cluster, and launching of MuTect jobs.
'''

import os
import re
from hashlib import md5

from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.pipeline.bwa_runner import ClusterJobManager, genome_fasta_path
from osqpipe.models import MergedAlnfile
from osqutil.config import Config

def cluster_fname(localfile):
  '''
  Generate a unique filename to be used on the cluster, without
  unnecessarily divulging local file paths.
  '''
  pathbits = os.path.split(localfile)
  hasher   = md5()
  hasher.update(pathbits[0])
  clusterfile = "%d_%s_%s" % (os.getpid(), hasher.hexdigest(), pathbits[1])
  return clusterfile

class MutectManager(ClusterJobManager):
  '''
  Class used to manage the transfer of input files to the cluster,
  and the submission of jobs to run MuTect, tranfer output files
  back to the local host, and clean up input and output files on the cluster.
  '''
  def __init__(self, genome, memsize=4, coverage=False,
               inprefixes=('IR_BQSR_ear_exome_','IR_BQSR_HCC_nodule_exome_'),
               **kwargs):

    super(MutectManager, self).__init__(memsize=memsize, **kwargs)
    self.config     = Config()
    self.genome     = genome
    self.coverage   = coverage
    self.inprefixes = inprefixes

  def _extract_stem(self, fname):
    '''
    Strip of filename suffix and any supplied prefix to get to some
    kind of core identifier.
    '''
    fname = os.path.split(fname)[-1]
    stem  = os.path.splitext(fname)[0]
    if self.inprefixes is not None:
      for pref in self.inprefixes:
        match = re.match(pref, fname)
        if match:
          return re.sub("^%s" % pref, '', stem)
    return stem

  def output_prefix(self, tumor, normal):
    '''
    Generate an output prefix which includes core IDs for both tumor and normal.
    '''
    return "MuTect_%s_%s" % (self._extract_stem(tumor), self._extract_stem(normal))

  def launch(self, tumor, normal):
    '''
    Main method used to initiate the MuTect cluster jobs.
    '''
    LOGGER.info("Copying files to cluster server.")

    # Note: .bai files need transfer also
    tumor_ix     = "%s.bai" % os.path.splitext(tumor)[0]
    normal_ix    = "%s.bai" % os.path.splitext(normal)[0]
    cl_tumor     = cluster_fname(tumor)
    cl_normal    = cluster_fname(normal)
    cl_tumor_ix  = cluster_fname(tumor_ix)
    cl_normal_ix = cluster_fname(normal_ix)
    self.submitter.remote_copy_files(filenames=(tumor,       normal,
                                                tumor_ix,    normal_ix),
                                     destnames=(cl_tumor,    cl_normal,
                                                cl_tumor_ix, cl_normal_ix))
    outprefix    = self.output_prefix(tumor, normal)
    cl_outprefix = cluster_fname(outprefix)

    # Determine the correct genome fasta path from self.genome
    # database Genome object.
    genome_fasta = genome_fasta_path(self.genome, self.config.clustergenomedir)

    LOGGER.info("Submitting MuTect SNV calling job.")
    cmd = (("mutect --analysis_type MuTect --reference_sequence %s"
            + " --input_file:normal %s --input_file:tumor %s"
            + " --out %s.out --vcf %s.vcf")
           % (genome_fasta, cl_normal, cl_tumor, cl_outprefix, cl_outprefix))

    # Optional output of wig file. This is big.
    if self.coverage:
      cmd += " --coverage_file %s_coverage.wig" % outprefix

    cl_queue = self.config.clusterqueue
    mutect_job = self.submit_command(cmd=cmd,
                                     mem=self.memsize * 1024,
                                     queue=cl_queue,
                                     auto_requeue=False)

    LOGGER.info("Submitting cleanup job.")
    cmd = "rm %s %s %s %s" % (cl_tumor, cl_normal, cl_tumor_ix, cl_normal_ix)
    self.submit_command(cmd=cmd,
                        queue=cl_queue,
                        depend_jobs=[mutect_job],
                        auto_requeue=False)

    # Transfer back to localhost and clean up outputs. Don't forget
    # the .vcf.idx.
    outpatterns = ("%s.out", "%s.vcf", "%s.vcf.idx")
    clustfiles  = [ pattern % cl_outprefix for pattern in outpatterns ]
    localfiles  = [ pattern % outprefix    for pattern in outpatterns ]
    LOGGER.info("Submitting output transfer jobs.")
    for num in range(len(outpatterns)):
      sshcmd = self.return_file_to_localhost(clustfiles[num], localfiles[num],
                                             execute=False)
      cmd = '%s && rm %s' % (sshcmd, clustfiles[num])
      self.submit_command(cmd=cmd,
                          queue=cl_queue,
                          depend_jobs=[mutect_job],
                          auto_requeue=False)

if __name__ == '__main__':

  from argparse import ArgumentParser
  from osqpipe.models import Genome

  PARSER = ArgumentParser(\
    description="Launch MuTect SNV calling on the cluster.")

  PARSER.add_argument('--tumour', '-t', dest='tumor', type=str, required=True,
                      help="The name of the tumour tissue bam file.")

  PARSER.add_argument('--normal', '-n', dest='normal', type=str, required=True,
                      help="The name of the normal tissue bam file.")

  PARSER.add_argument('--genome', '-g', dest='genome', type=str, required=True,
                      help="The code for the reference genome.")

  ARGS = PARSER.parse_args()

  GENOBJ = Genome.objects.get(code=ARGS.genome)

  MUTMAN = MutectManager(genome=GENOBJ)

  MUTMAN.launch(ARGS.tumor, ARGS.normal)
