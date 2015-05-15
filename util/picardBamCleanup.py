#!/usr/bin/env python

'''
Script to run picard CleanSam, AddOrReplaceReadGroups and
FixMateInformation on bam files stored in the repository, replacing
the old bam file with the cleaned-up version (the file name remains
the same; stored MD5 sums are updated)
'''

import os
from shutil import move

from osqpipe.pipeline.setup_logs import configure_logging
from osqpipe.pipeline.config import Config
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

from osqpipe.models import Alnfile
from osqpipe.pipeline.utilities import BamPostProcessor, call_subprocess, \
    set_file_permissions, checksum_file
from osqpipe.pipeline.samtools import count_bam_reads

CONFIG = Config()

def run_picard(libcode, facility, lanenum=None):

  if lanenum is None:
    bams = Alnfile.objects.filter(alignment__lane__library__code=libcode,
                                  alignment__lane__facility__code=facility,
                                  filetype__code='bam')
  else:
    bams = Alnfile.objects.filter(alignment__lane__library__code=libcode,
                                  alignment__lane__facility__code=facility,
                                  lane__lanenum=lanenum,
                                  filetype__code='bam')

  if len(bams) == 0:
    raise StandardError("Unable to find matching bam file in the repository.")

  for bam in bams:
    oldsum   = checksum_file(bam.repository_file_path)
    if oldsum != bam.checksum:
      raise ValueError(("MD5 checksum of bam file on disk (%s) does not agree"
                        + " with stored repository value (%s): %s")
                       % (oldsum, bam.checksum, bam.filename))
    newbam   = bam.repository_file_path + '.cleaned'
    postproc = BamPostProcessor(input_fn=bam.repository_file_path,
                                output_fn=newbam)

    # Run CleanSam
    call_subprocess(postproc.clean_sam(), path=CONFIG.hostpath)

    # Run AddOrReplaceReadGroups
    call_subprocess(postproc.add_or_replace_read_groups(), path=CONFIG.hostpath)
    os.unlink(postproc.cleaned_fn)

    # Run FixMateInformation
    call_subprocess(postproc.fix_mate_information(), path=CONFIG.hostpath)
    os.unlink(postproc.rgadded_fn)

    # Quick sanity check on the output
    newcount = count_bam_reads(newbam)

    # FIXME total_reads should be total reads in bam, not in fastq.
    if newcount != bam.alignment.total_reads:
      raise ValueError(("Read count in cleaned bam file (%d) does not agree"
                        + " with total_reads in repository (%d): %s")
                       % (newcount, bam.alignment.total_reads, newbam))

    # Clean up and replace the old bam file with the new one.
    set_file_permissions(CONFIG.group, newbam)
    checksum = checksum_file(newbam)
    os.unlink(bam.repository_file_path)
    move(newbam, bam.repository_file_path)
    bam.checksum = checksum
    bam.save()

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(
    description='Trim small RNA reads to remove adapter sequences and count exact read matches.')

  PARSER.add_argument('-l', '--library', dest='libcode', type=str, required=True,
                      help='The code of the library to process.')

  PARSER.add_argument('-f', '--facility', dest='facility', type=str, default='CRI',
                      help='The code of the sequencing facility (default=CRI)')

  PARSER.add_argument('-n', '--lanenum', dest='lanenum', type=int, required=False,
                      help='The (optional) lane number to which processing should be restricted.')

  ARGS = PARSER.parse_args()

  run_picard(ARGS.libcode, ARGS.facility, ARGS.lanenum)
    
