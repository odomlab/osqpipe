#!/usr/bin/env python
#
# $Id$

'''Given the name of a bam file generated by e.g. cs_processFile.py,
generate bed, wiggle and bedgraph files and store them all in the
repository.'''

import sys
import os
import re
from os.path import splitext
import logging
from pipes import quote
from subprocess import CalledProcessError
from shutil import move
from tempfile import NamedTemporaryFile

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

from osqpipe.pipeline.utilities import call_subprocess, \
    set_file_permissions, rezip_file, parse_repository_filename
from osqpipe.pipeline.config import Config
from osqpipe.models import Filetype, Genome, Library, \
    Status, Lane

from osqpipe.pipeline.samtools import BamToBedConverter
from osqpipe.pipeline.alignment import AlignmentHandler

###############################################################################
#
# Software dependencies:
# - bam2bed
# - makeWiggle
#
###############################################################################

###############################################################################

class AlignProcessingManager(object):

  '''
  Simple class encapsulating the creation and registration of bed,
  bedgraph and wiggle files.
  '''

  __slots__ = ('debug', 'conf')

  def __init__(self, debug=False):
    self.debug = debug
    if self.debug:
      LOGGER.setLevel(logging.DEBUG)
    else:
      LOGGER.setLevel(logging.INFO)
    self.conf = Config()

  def _reallocate_reads(self, in_fn):
    '''
    Run the reallocateReads script, overwriting the input file with a
    new bam file in which non-unique reads are reallocated according
    to the distribution of unique reads. A pretty horrible hack,
    statistically speaking.
    '''
    # Re-distribute non-unique reads.
    tmpfile = os.path.join(self.conf.tmpdir, "%s_cs_processAlignmentBwa.tmp" % in_fn)
    LOGGER.info("Re-allocating non-unique reads to a temporary file %s.", tmpfile)
    cmd = (self.conf.read_reallocator, in_fn, tmpfile)
    LOGGER.debug(cmd)
    call_subprocess(cmd, path=self.conf.hostpath)

    # Re-sort output bam
    LOGGER.info("Sorting temporary file %s back to %s.", tmpfile, in_fn)
    cmd = (self.conf.read_sorter, 'sort', '-m', self.conf.meminbytes, tmpfile, splitext(in_fn)[0])
    LOGGER.debug(cmd)
    call_subprocess(cmd, path=self.conf.hostpath)

    # Remove tmp file
    LOGGER.info("Removing temporary file %s.", tmpfile)
    os.unlink(tmpfile)

  def get_genome_size_file(self, genome):
    '''
    Retrieve the filename containing chromosome lengths for a given
    genome. Returns two values: the name of the filename, and whether
    that filename should be treated as a temporary file, i.e. to be
    deleted once done with. Such deletion is the responsibility of the
    calling code.
    '''
    fnchrlen = os.path.join(self.conf.genomesizedir, genome + ".fa.length")
    if not os.path.exists(fnchrlen):
      tmpfile = NamedTemporaryFile(delete=False, dir=self.conf.tmpdir)
      cmd = "%s %s > %s" % ('fetchChromSizes',
                            genome,
                            tmpfile.name)
  
      # Note - assumes we're running on our primary host. FIXME?
      call_subprocess(cmd, shell=True, path=self.conf.hostpath)
      tmpfile.close()
      try:
        LOGGER.info("Storing new chromosome sizes file as %s", fnchrlen)
        move(tmpfile.name, fnchrlen)
      except Exception, err:
        LOGGER.warning("Attempt to store chromosome sizes file"
                       + " as %s failed: %s", fnchrlen, err)
        return(tmpfile.name, True)
    else:
      LOGGER.debug("Retrieved cached chromosome sizes file %s", fnchrlen)
    return (fnchrlen, False)

  def check_bam_vs_lane_fastq(self, bam, relaxed=False):
    '''
    Quick check that the number of reads in the bam file being saved
    is identical to the number of (passed PF) reads in the input fastq
    file.
    '''
    (code, facility, lanenum, _pipeline) = parse_repository_filename(bam)
    try:
      lane = Lane.objects.get(library__code=code,
                              facility__code=facility,
                              lanenum=lanenum)
    except Lane.DoesNotExist, err:
      LOGGER.error("Unexpected lane in filename, not found in repository.")
      sys.exit("Unable to find lane in repository")

    LOGGER.info("Checking number of reads in bam file %s", bam)
    cmd  = (self.conf.read_sorter, 'flagstat', bam)
    pout = call_subprocess(cmd, path=self.conf.hostpath)
    numreads = int(pout.readline().split()[0])
    expected = lane.passedpf
    if lane.paired:
      expected *= 2
    if numreads != expected:
      message = ("Number of reads in bam file is differs from that in "
                 + "fastq file: %d (bam) vs %d (fastq)")
      if relaxed:
        LOGGER.warning(message, numreads, expected)
      else:
        raise ValueError(message % (numreads, expected))
    
  def run(self, in_fn, genome, reallocate, relaxed=False):

    '''Given a bam file name, run the file conversions and store all
    files in the repository (input bam file included).'''

    bedtype = Filetype.objects.get(code='bed')

    base = splitext(in_fn)[0]

    try:
      lib = Library.objects.search_by_filename(in_fn)
    except Library.DoesNotExist, err:
      LOGGER.error("Unable to determine library from filename %s" % (in_fn,))
      sys.exit("Unexpected file naming convention.")
      
    if genome is None:
      genome = lib.genome.code

    # The actual database loading code is currently set to raise an
    # error in such cases. An earlier warning would be nice, though.
    if not re.search(genome, in_fn):
      LOGGER.warning("Filename does not match expected genome (%s) and so database loading may fail." % genome)

    self.check_bam_vs_lane_fastq(in_fn, relaxed)

    # Set aligner, later passed as argument to AlignmentHandler.
    aligner = self.conf.aligner
    params  = ''
  
    # In the case of PolIII/TFIIIC libraries, assume that the bam file
    # was created by keeping many non-unique reads. Note: we assume
    # here that re-running the reallocateReads script over the output
    # file would have little to no effect, and so we do not wrap this
    # within our database transaction.
    if reallocate or lib.libtype.name == "ChIP-Seq" \
          and lib.factor is not None \
          and lib.factor.name in self.conf.reallocation_factors:

      self._reallocate_reads(in_fn)

      # Set the aligner value for later.
      aligner = [ self.conf.aligner, self.conf.read_reallocator, self.conf.read_sorter ]

      # This would have been previously set in BwaClusterJobSubmitter FIXME DRY.
      params  = [ '--n_occ %s' % self.conf.nonuniquereads, '', '' ]

    (chrom_sizes, chr_istmp) = self.get_genome_size_file(genome)

    # make bed file(s)
    bed_fn = base + bedtype.suffix
    beds  = []
    bamToBed = BamToBedConverter(tc1         = genome == 'tc1',
                                 chrom_sizes = chrom_sizes)
    for outBed in bamToBed.convert(in_fn, bed_fn):
      beds.append(outBed)

    # make wiggle files
    wigs = []
    for bed_fn in beds:
      cmd = ('makeWiggle', '-t', '3', bed_fn, splitext(bed_fn)[0])
      LOGGER.debug(" ".join(cmd))
      pout = call_subprocess(cmd, path=self.conf.hostpath)
      for line in pout:
        wigs.append(line.strip())
      pout.close()

    # make bedgraph file(s)
    LOGGER.info("Creating bedgraph files...")
    bedgraphs = []
    for bed_fn in beds:
      cmd = ('makeWiggle', '-t', '3', '-B', '-1', bed_fn, splitext(bed_fn)[0])
      LOGGER.debug(cmd)
      try:  # Occasionally fails for poorer-quality genomes.
        pout = call_subprocess(cmd, path=self.conf.hostpath)
        for line in pout:
          bedgraphs.append(line.strip())
        pout.close()
      except CalledProcessError, err:
        LOGGER.warn("Unable to create bedGraph file: %s", err)

    # make bigWig file(s)
    LOGGER.info("Creating bigWig files...")
    bigwigs = []
    for bgr_fn in bedgraphs:
      bwfile = splitext(bgr_fn)[0] + '.bw'
      cmd = ('bedGraphToBigWig', bgr_fn, chrom_sizes, bwfile)
      LOGGER.debug(cmd)
      try:  # This can fail, e.g. for very small input files.
        call_subprocess(cmd, path=self.conf.hostpath)
        bigwigs.append(bwfile)
      except CalledProcessError, err:
        LOGGER.warn("Unable to create bigWig file: %s", err)
        
    # We're now done with the chrom_sizes file.
    if chr_istmp:
      os.unlink(chrom_sizes)

    # Set group ownership and permissions appropriately
    grp = self.conf.group
    set_file_permissions(grp, in_fn)
    for bed in beds:
      set_file_permissions(grp, bed)
    for wig in wigs:
      set_file_permissions(grp, wig)
    for bgr in bedgraphs:
      set_file_permissions(grp, bgr)
    for bwig in bigwigs:
      set_file_permissions(grp, bwig)
  
    # compress bed file(s)
    bedgz = []
    for bed in beds:
      gzname = rezip_file(bed)
      bedgz.append(gzname)

    # compress wiggle files
    wigsgz = []
    for wig in wigs:
      gzname = rezip_file(wig)
      wigsgz.append(gzname)

    # compress bedgraph files
    bgrgz = []
    for bgr in bedgraphs:
      gzname = rezip_file(bgr)
      bgrgz.append(gzname)

    # Don't compress bigwig files

    # Add alignment to database.  The AlignmentHandler class has its
    # own transaction, however we want to make sure that the lane
    # status is updated as well. This is best done by passing the
    # desired final status into the transaction itself. We don't
    # simply wrap this in a nested transaction because the
    # AlignmentHandler computes file checksums (outside of its
    # transaction) and that takes time.
    hnd     = AlignmentHandler(prog=aligner, params=params, genome=genome)
    fstatus = Status.objects.get(code='complete', authority=None)
    lane    = hnd.add(bedgz + bgrgz + wigsgz + bigwigs + [in_fn],
                      final_status=fstatus)

    # Occasionally we process a bam file we're confident is completely
    # transferred but where we've forgotten to create a bam.done file. Such
    # omissions shouldn't crash the script.
    done_fn = in_fn + ".done"
    if os.path.exists(done_fn):
      LOGGER.debug("rm -f %s" % (done_fn,))
      os.unlink(done_fn)
    else:
      LOGGER.warn("File not found for scheduled deletion: %s", done_fn)

#######################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Various postprocessing steps for alignments.')

  PARSER.add_argument('file', metavar='<filename.bam>', type=str,
                      help='The bam file name (required).')

  PARSER.add_argument('-d', '--debug', dest='debug', action='store_true',
                      help='Turn on debugging mode.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=False,
                      help='The genome used for the alignment. Use only if different from one provided in library.')

  PARSER.add_argument('-r', '--reallocate', dest='reallocate', action='store_true', required=False,
                      help='Re-allocates non-unique reads.')

  PARSER.add_argument('--relaxed', dest='relaxed', action='store_true', required=False,
                      help='Ignore validation errors (e.g., mismatches'
                      + ' between numbers of reads in bam and fastq files.')

  ARGS = PARSER.parse_args()

  APM = AlignProcessingManager(debug=ARGS.debug)

  APM.run(ARGS.file, ARGS.genome, ARGS.reallocate, ARGS.relaxed)
