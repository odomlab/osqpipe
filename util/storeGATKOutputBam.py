#!/usr/bin/env python

'''
Script to take the final outputs of runGATKPreprocessing.py and store
them in our repository as MergedAlignment/MergedAlnfile pairs.
'''

import os
import re

from shutil import move
from logging import INFO
from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)
CONFIG = Config()

# New in Django 1.7 and above.
import django
django.setup()

from pysam import AlignmentFile
from django.db import transaction
from osqpipe.models import MergedAlignment, MergedAlnfile, Filetype
from osqpipe.pipeline.utilities import checksum_file, set_file_permissions
from osqpipe.pipeline.config import Config
from osqpipe.pipeline.gatk import retrieve_readgroup_alignment, check_bam_readcount

@transaction.atomic
def load_merged_bam(bam, genome=None):

  LOGGER.info("Storing bam file %s in repository.", bam)

  bamtype = Filetype.objects.get(code='bam')
  
  with AlignmentFile(filename=bam) as bamhandle:
    rgroups = bamhandle.header.get('RG', [])

  samples = list(set([ rg.get('SM') for rg in rgroups ]))
  if len(samples) > 1:
    raise StandardError("More than one sample represented in input.")
  elif len(samples) == 0:
    LOGGER.warning("No samples specified in bam file.")

  # Slightly convoluted multiple query (as opposed to query__in) so
  # we can be sure we're identifying all the target Alignments.
  alns = [ retrieve_readgroup_alignment(rg, genome) for rg in rgroups ]
  alns = list(set(alns))

  LOGGER.info("Linking MergedAlignment to %d source Alignments.", len(alns))
  maln = MergedAlignment.objects.create()
  for aln in alns:
    maln.alignments.add(aln)
  maln.full_clean() # Raise ValidationError if the MergedAlignment contains inconsistencies.

  LOGGER.info("Calculating bam file MD5 checksum.")
  chksum = checksum_file(bam, unzip=False)

  LOGGER.info("Counting reads in bam file.")
  check_bam_readcount(bam, maln)

  malnfile = MergedAlnfile.objects.create(alignment=maln,
                                          filename=bam,
                                          filetype=bamtype,
                                          checksum=chksum)

  LOGGER.info("Moving file into repository.")
  destname = malnfile.repository_file_path
  move(bam, destname)
  set_file_permissions(CONFIG.group, destname)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(\
    description="Store GATK pipeline output files in the repository.")

  PARSER.add_argument('bams', metavar='<bam files>', type=str, nargs='+',
                      help='The names of the bam files to store in the repository.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=False,
                      help='The genome used to filter the alignments against which'
                      + ' the merged files may be linked. Use in cases where there'
                      + ' is ambiguity; otherwise the entire pool of alignments for'
                      + ' the lanes of interest will be available for linking.')

  ARGS = PARSER.parse_args()

  for bam in ARGS.bams:
    load_merged_bam(bam, ARGS.genome)

    # Remove bam.done file if present.
    donefile = "%s.done" % bam
    if os.path.exists(donefile):
      os.unlink(donefile)
