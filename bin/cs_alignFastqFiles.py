#!/usr/bin/env python
#
# $Id$

'''Script to run the alignment of a given fastq file or files against
a genome registered in the repository.'''

from osqpipe.pipeline.fastq_aligner import FastqBwaAligner
from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Align one or a pair of Fastq files against a specified genome.')

  PARSER.add_argument('files', metavar='<filenames>', type=str, nargs='+',
                      help='The name of the file or files to process.\n'
                      + 'The files may be a single or a pair of fastq files.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=True,
                      help='The genome used in the alignment.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Turn on test mode.')

  PARSER.add_argument('--nocleanup', dest='nocleanup', action='store_true',
                      help='Keeps all temporary files.')

  PARSER.add_argument('--n_occ', dest='nocc', type=str,
                      help='Specifies number of non-unique read occurrences to keep in bam file. The bwa default is 3.')

  ARGS = PARSER.parse_args()

  BWA = FastqBwaAligner(testMode=ARGS.testMode)
  
  BWA.align_standalone(filepaths=ARGS.files,
                       genome  = ARGS.genome,
                       nocleanup = ARGS.nocleanup,
                       nocc = ARGS.nocc)
  


