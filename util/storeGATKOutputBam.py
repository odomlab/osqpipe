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
Script to take the final outputs of runGATKPreprocessing.py and store
them in our repository as MergedAlignment/MergedAlnfile pairs.
'''

import os
import time

from shutil import move
from logging import INFO
from osqutil.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from pysam import AlignmentFile
from django.db import transaction
from osqpipe.models import MergedAlignment, MergedAlnfile, Filetype, \
    Alignment, ArchiveLocation, Program, DataProvenance
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
      if read.is_supplementary:
        continue
      rgid = read.get_tag('RG')
      total[rgid] = total.setdefault(rgid, 0) + 1
      if not read.is_unmapped:
        mapped[rgid] = mapped.setdefault(rgid, 0) + 1
        if read.mapq > 0:
          munique[rgid] = munique.setdefault(rgid, 0) + 1

  return (total, mapped, munique)

@transaction.atomic
def _store_bam_within_transaction(bam, rgroups, readcountdicts,
                                  genome=None, bamfilter=False, autoaln=False,
                                  alignprog=None, alignparams=None,
                                  archloc=None):
  '''
  Function to make all the database changes within a single short
  transaction.
  '''
  bamtype = Filetype.objects.get(code='bam')

  # Slightly convoluted multiple query (as opposed to query__in) so we
  # can be sure we're identifying all the target Alignments. We also
  # auto-create Alignments here if necessary.
  alns = []
  for rgp in rgroups:
    try:
      alns += [ retrieve_readgroup_alignment(rgp, genome, bamfilter) ]
    except Alignment.DoesNotExist, err:
      if autoaln:
        newaln = autocreate_alignment(rgp, genome,
                                      [ rgcount.get(rgp.get('ID'), 0)
                                        for rgcount in readcountdicts ])
        if alignprog is not None:
          DataProvenance.objects.create(program      = alignprog,
                                        parameters   = alignparams,
                                        rank_index   = 0,
                                        data_process = newaln)          
        alns += [ newaln ]
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
  check_bam_readcount(bam, maln, readcountdicts)

  malnfile = MergedAlnfile.objects.create(alignment=maln,
                                          filename=bam,
                                          filetype=bamtype,
                                          checksum=chksum)

  if archloc is not None:
    malnfile.archive = ArchiveLocation.objects.get(name=archloc)
    malnfile.archive_date = time.strftime('%Y-%m-%d')
    malnfile.save()

  LOGGER.info("Moving file into repository.")
  destname = malnfile.repository_file_path
  move(bam, destname)
  set_file_permissions(CONFIG.group, destname)

def load_merged_bam(bam, genome=None, bamfilter=False, autoaln=False,
                    aligner=None, alignvers=None, alignparams=None,
                    archloc=None):
  '''
  Insert the specified merged bam file into the repository, linking
  against per-lane Alignments as appropriate.
  '''
  if archloc is None:
    LOGGER.info("Storing merged bam file %s in repository...", bam)
  else:
    LOGGER.info("Storing merged bam file %s in archive %s...", bam, archloc)

  alignprog = None
  if aligner is not None:
    if alignvers is None:
      raise ValueError("Aligner version must be set if specifying aligner program")
    try:
      alignprog = Program.objects.get(program=aligner,
                                      version=alignvers,
                                      current=True)
    except Program.DoesNotExist, _err:
      raise StandardError("Unable to find current program in database: %s %s"
                          % (aligner, alignvers))

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
  readcountdicts = None
  if autoaln:
    LOGGER.info("Counting reads per readgroup for %s...", bam)
    readcountdicts = count_readgroup_reads(bam)

  _store_bam_within_transaction(bam, rgroups, readcountdicts,
                                genome, bamfilter, autoaln,
                                alignprog, alignparams,
                                archloc)

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

  PARSER.add_argument('--aligner', dest='aligner', type=str, required=False,
                      help='The aligner with which to annotate any automatically'
                      + ' created Alignment objects.')

  PARSER.add_argument('--aligner-version', dest='alignvers', type=str, required=False,
                      help='The aligner version with which to annotate any automatically'
                      + ' created Alignment objects.')

  PARSER.add_argument('--aligner-params', dest='alignparams', type=str, required=False,
                      help='The aligner parameters with which to annotate any automatically'
                      + ' created Alignment objects.')

  PARSER.add_argument('--archive', dest='archloc', type=str, required=False,
                      help='The archive location destination in which to store'
                      + ' files (defaults to the configured location).')

  ARGS = PARSER.parse_args()

  for filename in ARGS.bams:
    load_merged_bam(filename,
                    ARGS.genome,
                    ARGS.bamfilt,
                    ARGS.autoaln,
                    ARGS.aligner,
                    ARGS.alignvers,
                    ARGS.alignparams,
                    ARGS.archloc)

    # Remove bam.done file if present.
    donefile = "%s.done" % filename
    if os.path.exists(donefile):
      os.unlink(donefile)
