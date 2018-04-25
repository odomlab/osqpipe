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
Script to run the Kraken reaper and tally tools on the files for a
given library. This may be used to update the files we keep for older
smallRNA libraries. Note that you will need to insert the new files
into the repository manually.
'''

import os

from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane
from osqpipe.pipeline.file_processor import MiRFastqFileProc

def process_smallrna_fastq(fastq, lane):

  proc = MiRFastqFileProc(fname=fastq.repository_file_path,
                          fname2=None,
                          paired=lane.paired,
                          facility=lane.facility,
                          flowcell=lane.flowcell,
                          flowlane=lane.flowlane,
                          libcode=lane.library.code)

  screened  = proc.trim_linkers(fastq.repository_file_path)
  proc.tempfiles.append(screened)
  clustered = proc.cluster_exact_matches(screened,
                                         os.path.splitext(fastq.filename)[0])
  proc.clean_up()

def process_smallrna_lane(lane):

  LOGGER.info("Processing lane %s" % lane)

  files = lane.lanefile_set.filter(filetype__code='fq')

  if len(files) not in (1,2):
    raise ValueError("Unexpected number of fastq files returned by the repository.")

  for fastq in files:
    process_smallrna_fastq(lane=lane,
                           fastq=fastq)

def process_smallrna_library(libcode, lanenum=None):

  if lanenum is None:
    lanes = Lane.objects.filter(library__code=libcode)
  else:
    lanes = Lane.objects.filter(library__code=libcode, lane__lanenum=lanenum)

  if len(lanes) == 0:
    raise StandardError("Unable to find library (and lane) in the repository.")

  for lane in lanes:
    process_smallrna_lane(lane)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(
    description='Trim small RNA reads to remove adapter sequences and count exact read matches.')

  PARSER.add_argument('-l', '--library', dest='libcode', type=str, required=True,
                      help='The name of the library to process.')

  PARSER.add_argument('-n', '--lanenum', dest='lanenum', type=int, required=False,
                      help='The (optional) lane number to which processing should be restricted.')

  ARGS = PARSER.parse_args()

  process_smallrna_library(ARGS.libcode, ARGS.lanenum)
    
