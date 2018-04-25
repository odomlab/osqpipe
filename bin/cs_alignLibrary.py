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

'''Script to rerun the alignment of a given repository library's fastq
file against a genome registered in the repository.'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.pipeline.fastq_aligner import FastqBwaAligner, FastqTophatAligner, FastqStarAligner
from osqpipe.models import Library

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Align a given library against a specified genome.')

  PARSER.add_argument('-l', '--library', dest='library', type=str, required=True,
                      help='The code for the library to align.')

  PARSER.add_argument('-f', '--facility', dest='facility', type=str, required=False,
                      help='The facility in which the lane was sequenced. If neither this'
                      + ' nor lanenum are specified, all lanes of the library are aligned.')

  PARSER.add_argument('-n', '--lanenum', dest='lanenum', type=str, required=False,
                      help='The lane number of the library. If not specified, all'
                      + ' lanes of the library are aligned.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=True,
                      help='The genome used in the alignment.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Turn on test mode.')

  PARSER.add_argument('--nocleanup', dest='nocleanup', action='store_true',
                      help='Keeps all temporary files.')

  PARSER.add_argument('--n_occ', dest='nocc', type=str,
                      help='Specifies number of non-unique read occurrences to keep'
                      + ' in bam file. The bwa default is 3.')

  PARSER.add_argument('--aligner', type=str, dest='aligner', choices=('bwa', 'tophat', 'star'), default=None, 
                      help='Tophat is the default aligner for rnaseq data while bwa is default for all other library types.')
  
  PARSER.add_argument('--algorithm', type=str, dest='algorithm', choices=('aln', 'mem'),
                      help='The bwa algorithm to use (aln or mem). The default behaviour'
                      + ' is to pick the algorithm based on the read length in the fastq files.')

  ARGS = PARSER.parse_args()

  # Using a subclass to define the external program will allow us to
  # tailor this to other aligners in future.
  library = Library.objects.get(code=ARGS.library)
  if library.libtype.code == 'rnaseq':
    if ARGS.aligner in (None, 'tophat'):                 # The default.
      BWA = FastqTophatAligner(test_mode=ARGS.testMode,
                               samplename=library.sample.name)
    elif ARGS.aligner == 'star':
      BWA = FastqStarAligner(test_mode=ARGS.testMode,
                             samplename=library.sample.name)
    else:
      raise StandardError("Unrecognised aligner for %s: %s" % (library.libtype.name, ARGS.aligner))
  elif ARGS.aligner in (None, 'bwa'):
    BWA = FastqBwaAligner(test_mode=ARGS.testMode,
                          samplename=library.sample.name,
                          bwa_algorithm=ARGS.algorithm)
  else:
    raise StandardError("Unrecognised aligner for %s: %s" % (library.libtype.name, ARGS.aligner))

  BWA.align(library = ARGS.library,
            facility = ARGS.facility,
            lanenum = ARGS.lanenum,
            genome  = ARGS.genome,
            nocleanup = ARGS.nocleanup,
            nocc = ARGS.nocc)
