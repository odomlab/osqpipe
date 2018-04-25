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

'''Script to run the alignment of a given fastq file or files against
a genome registered in the repository.'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.pipeline.fastq_aligner import FastqBwaAligner, FastqTophatAligner, FastqStarAligner

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Align one or a pair of Fastq files against a specified genome.')

  PARSER.add_argument('files', metavar='<filenames>', type=str, nargs='+',
                      help='The name of the file or files to process.\n'
                      + 'The files may be a single or a pair of fastq files.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=True,
                      help='The genome used in the alignment.')

  PARSER.add_argument('-s', '--sample', dest='sample', type=str, required=False,
                      help='The sample name with which to tag the output bam read groups.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Turn on test mode.')

  PARSER.add_argument('--nocleanup', dest='nocleanup', action='store_true',
                      help='Keeps all temporary files.')

  PARSER.add_argument('--n_occ', dest='nocc', type=str,
                      help='Specifies number of non-unique read occurrences to keep'
                      + ' in bam file. The bwa default is 3.')

  PARSER.add_argument('-a', '--aligner', dest='aligner', type=str,
                      choices=('bwa', 'tophat', 'star'), default='bwa',
                      help='The aligner program to use.')

  PARSER.add_argument('--algorithm', type=str, dest='algorithm', choices=('aln', 'mem'),
                      help='The bwa algorithm to use (aln or mem). The default behaviour'
                      + ' is to pick the algorithm based on the read length in the fastq files.')

  PARSER.add_argument('--rcp', type=str, dest='rcp',
                      help='Remote file copy (rcp) target.')

  PARSER.add_argument('--lcp', type=str, dest='lcp', default=None,
                      help='Local file copy (lcp) target.')
  
  PARSER.add_argument('--no-split', dest='nosplit', action='store_true',
                      help='Do not split input fastq for distributed parallel alignment.', default=False)

  PARSER.add_argument('--fileshost', dest='fileshost', type=str,
                      help='Host where the files should be downloaded from.')


  ARGS = PARSER.parse_args()

  if ARGS.aligner == 'bwa':
    BWA = FastqBwaAligner(test_mode=ARGS.testMode,
                          samplename=ARGS.sample,
                          bwa_algorithm=ARGS.algorithm)
  elif ARGS.aligner == 'tophat':
    BWA = FastqTophatAligner(test_mode=ARGS.testMode,
                             samplename=ARGS.sample)
  elif ARGS.aligner == 'star':
    BWA = FastqStarAligner(test_mode=ARGS.testMode,
                             samplename=ARGS.sample)
  else:
    raise ValueError("Unrecognised aligner requested: %s" % ARGS.aligner)
  
  BWA.align_standalone(filepaths=ARGS.files,
                       genome  = ARGS.genome,
                       nocleanup = ARGS.nocleanup,
                       nocc = ARGS.nocc,
                       nosplit=ARGS.nosplit,
                       rcp=ARGS.rcp,
                       lcp=ARGS.lcp,
                       fileshost=ARGS.fileshost)
