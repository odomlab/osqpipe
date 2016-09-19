#!/usr/bin/env python

'''
Script to take the final outputs of runGATKPreprocessing.py and store
them in our repository as MergedAlignment/MergedAlnfile pairs.
'''

import os

from shutil import move
from logging import INFO
from osqutil.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from pysam import AlignmentFile
from django.db import transaction
from osqpipe.models import MergedAlignment, MergedAlnfile, Filetype, Alignment
from osqutil.utilities import checksum_file, set_file_permissions
from osqutil.config import Config
from osqpipe.pipeline.gatk import retrieve_readgroup_alignment, \
    check_bam_readcount, autocreate_alignment
from osqpipe.pipeline.bampy import open_bamfile

CONFIG = Config()

def count_readgroup_reads(bam):
  '''
  Count the mapped and munique reads for each read group in a bam file.
  '''
  # This is liable to create a bam index file if one is not already
  # available.
  total   = dict()
  mapped  = dict()
  munique = dict()
  with open_bamfile(bam) as bamfile:
    for read in bamfile.fetch(until_eof=True):
      rgid = read.tags.get('RG')
      total[rgid] = total.setdefault(rgid, 0) + 1
      if not read.is_unmapped:
        mapped[rgid] = mapped.setdefault(rgid, 0) + 1
        if read.mapq > 0:
          munique[rgid] = munique.setdefault(rgid, 0) + 1

  return (total, mapped, munique)

@transaction.atomic
def load_merged_bam(bam, genome=None, bamfilter=False, autoaln=False):
  '''
  Insert the specified merged bam file into the repository, linking
  against per-lane Alignments as appropriate.
  '''
  LOGGER.info("Storing merged bam file %s in repository...", bam)

  bamtype = Filetype.objects.get(code='bam')

  # We use AlignmentFile rather than open_bamfile where possible,
  # because the latter will typically index the bam before doing
  # anything and that's not always desirable.
  with AlignmentFile(filename=bam) as bamhandle:
    rgroups = bamhandle.header.get('RG', [])

  samples = list(set([ rg.get('SM') for rg in rgroups ]))
  if len(samples) > 1:
    raise StandardError("More than one sample represented in input.")
  elif len(samples) == 0:
    LOGGER.warning("No samples specified in bam file.")

  # Count reads in readgroups for autoaln, if activated.
  readcountdict = None
  if autoaln:
    readcountdict = count_readgroup_reads(bam)

  # Slightly convoluted multiple query (as opposed to query__in) so we
  # can be sure we're identifying all the target Alignments. We also
  # auto-create Alignments here if necessary.
  alns = []
  for rgp in rgroups:
    try:
      alns += [ retrieve_readgroup_alignment(rgp, genome, bamfilter) ]
    except Alignment.DoesNotExist, err:
      if autoaln:
        autocreate_alignment(rgp, genome,
                             [ rgcount.get(rgp.get('ID'), 0)
                               for rgcount in readcountdict ])
      else:
        raise err

  alns = list(set(alns))

  LOGGER.info("Linking MergedAlignment to %d source Alignments.", len(alns))
  maln = MergedAlignment.objects.create()
  for aln in alns:
    maln.alignments.add(aln)

  # Raise ValidationError if the MergedAlignment contains inconsistencies.
  maln.full_clean()

  LOGGER.info("Calculating bam file MD5 checksum...")
  chksum = checksum_file(bam, unzip=False)

  LOGGER.info("Checking read count in bam file against lane records...")
  check_bam_readcount(bam, maln, readcountdict)

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
                      help='The names of the bam files to store in the'
                      + ' repository.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=False,
                      help='The genome used to filter the alignments against'
                      + ' which the merged files may be linked. Use in cases'
                      + ' where there is ambiguity; otherwise the entire pool'
                      + ' of alignments for the lanes of interest will be'
                      + ' available for linking.')

  PARSER.add_argument('--filter-alns-without-bams', dest='bamfilt',
                      action='store_true',
                      help='When linking the merged files to alignments, filter'
                      + ' out those alignments which do not include a bam file.'
                      + ' This is an occasionally useful workaround for some'
                      + ' edge cases.')

  PARSER.add_argument('--auto-create-alns', dest='autoaln', action='store_true',
                      help='Automatically create missing Alignment objects'
                      + ' using the data in the bam files.')

  ARGS = PARSER.parse_args()

  for filename in ARGS.bams:
    load_merged_bam(filename, ARGS.genome, ARGS.bamfilt, ARGS.autoaln)

    # Remove bam.done file if present.
    donefile = "%s.done" % filename
    if os.path.exists(donefile):
      os.unlink(donefile)
