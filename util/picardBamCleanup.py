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
Script to run picard CleanSam, AddOrReplaceReadGroups and
FixMateInformation on bam files stored in the repository, replacing
the old bam file with the cleaned-up version (the file name remains
the same; stored MD5 sums are updated)
'''

import os
from shutil import move

from osqutil.setup_logs import configure_logging
from osqutil.config import Config
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Alnfile
from osqutil.utilities import BamPostProcessor, call_subprocess, \
    set_file_permissions, checksum_file
from osqutil.samtools import count_bam_reads

from django.db import transaction

CONFIG = Config()

@transaction.atomic
def replace_repo_file(bam, newbam):

  bam = Alnfile.objects.get(id=bam.id) # Reload passed object within transaction.
  set_file_permissions(CONFIG.group, newbam)
  checksum = checksum_file(newbam, unzip=False)
  bam.checksum = checksum
  os.unlink(bam.repository_file_path)
  move(newbam, bam.repository_file_path)
  bam.save()

def run_picard(libcode, facility, lanenum=None, genome=None):

  bams = Alnfile.objects.filter(alignment__lane__library__code=libcode,
                                alignment__lane__facility__code=facility,
                                filetype__code='bam')
  if lanenum is not None:
    bams = bams.filter(alignment__lane__lanenum=lanenum)
  if genome is not None:
    bams = bams.filter(alignment__genome__code=genome)

  if len(bams) == 0:
    raise StandardError("Unable to find matching bam file in the repository.")

  for bam in bams:
    LOGGER.info("Confirming file checksum: %s", bam.filename)
    oldsum   = checksum_file(bam.repository_file_path, unzip=False)
    if oldsum != bam.checksum:
      raise ValueError(("MD5 checksum of bam file on disk (%s) does not agree"
                        + " with stored repository value (%s): %s")
                       % (oldsum, bam.checksum, bam.filename))
    newbam   = bam.repository_file_path + '.cleaned'
    postproc = BamPostProcessor(input_fn=bam.repository_file_path,
                                output_fn=newbam)

    # Run CleanSam
    LOGGER.info("Running CleanSam...")
    call_subprocess(postproc.clean_sam(), path=CONFIG.hostpath)

    # Run AddOrReplaceReadGroups
    LOGGER.info("Running AddOrReplaceReadGroups...")
    call_subprocess(postproc.add_or_replace_read_groups(), path=CONFIG.hostpath)
    os.unlink(postproc.cleaned_fn)

    # Run FixMateInformation
    LOGGER.info("Running FixMateInformation...")
    call_subprocess(postproc.fix_mate_information(), path=CONFIG.hostpath)
    os.unlink(postproc.rgadded_fn)

    # Quick sanity check on the output
    newcount = count_bam_reads(newbam)

    # FIXME total_reads should be total reads in bam, not in fastq.
    oldcount = bam.alignment.total_reads
    if bam.alignment.lane.paired:
      oldcount = oldcount * 2
    if newcount != oldcount:
      raise ValueError(("Read count in cleaned bam file (%d) does not agree"
                        + " with total_reads in repository (%d): %s")
                       % (newcount, oldcount, newbam))

    # Clean up and replace the old bam file with the new one.
    LOGGER.info("Replacing old bam file with new: %s", bam.repository_file_path)
    replace_repo_file(bam, newbam)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(
    description='Trim small RNA reads to remove adapter sequences and count exact read matches.')

  PARSER.add_argument('-l', '--library', dest='libcode', type=str, required=True,
                      help='The code of the library to process.')

  PARSER.add_argument('-f', '--facility', dest='facility', type=str, default='CRI',
                      help='The code of the sequencing facility (default=CRI).')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=False,
                      help='The genome code used for the alignments of interest.')

  PARSER.add_argument('-n', '--lanenum', dest='lanenum', type=int, required=False,
                      help='The (optional) lane number to which processing should be restricted.')

  ARGS = PARSER.parse_args()

  run_picard(ARGS.libcode, ARGS.facility, ARGS.lanenum, ARGS.genome)
    
