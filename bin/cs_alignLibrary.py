#!/usr/bin/env python
#
# $Id$

'''Script to rerun the alignment of a given repository library's fastq
file against a genome registered in the repository.'''

from osqpipe.pipeline.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

from osqpipe.pipeline.fastq_aligner import FastqBwaAligner, FastqTophatAligner
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
                      help='The lane number of the library. If not specified, all lanes of the library are aligned.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=True,
                      help='The genome used in the alignment.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Turn on test mode.')

  PARSER.add_argument('--nocleanup', dest='nocleanup', action='store_true',
                      help='Keeps all temporary files.')

  PARSER.add_argument('--n_occ', dest='nocc', type=str,
                      help='Specifies number of non-unique read occurrences to keep in bam file. The bwa default is 3.')

  ARGS = PARSER.parse_args()

  # Using a subclass to define the external program will allow us to
  # tailor this to other aligners in future.
  library = Library.objects.get(code=ARGS.library)
  if library.libtype.code == 'rnaseq':
    BWA = FastqTophatAligner(test_mode=ARGS.testMode, samplename=library.individual)
  else:
    BWA = FastqBwaAligner(test_mode=ARGS.testMode, samplename=library.individual)
  
  BWA.align(library = ARGS.library,
            facility = ARGS.facility,
            lanenum = ARGS.lanenum,
            genome  = ARGS.genome,
            nocleanup = ARGS.nocleanup,
            nocc = ARGS.nocc)
  


