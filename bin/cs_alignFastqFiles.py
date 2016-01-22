#!/usr/bin/env python
#
# $Id$

'''Script to run the alignment of a given fastq file or files against
a genome registered in the repository.'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.pipeline.fastq_aligner import FastqBwaAligner, FastqTophatAligner

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
                      choices=('bwa', 'tophat'), default='bwa',
                      help='The aligner program to use.')

  PARSER.add_argument('--algorithm', type=str, dest='algorithm', choices=('aln', 'mem'),
                      help='The bwa algorithm to use (aln or mem). The default behaviour'
                      + ' is to pick the algorithm based on the read length in the fastq files.')

  ARGS = PARSER.parse_args()

  if ARGS.aligner == 'bwa':
    BWA = FastqBwaAligner(test_mode=ARGS.testMode,
                          samplename=ARGS.sample,
                          bwa_algorithm=ARGS.algorithm)
  elif ARGS.aligner == 'tophat':
    BWA = FastqTophatAligner(test_mode=ARGS.testMode,
                             samplename=ARGS.sample)
  else:
    raise ValueError("Unrecognised aligner requested: %s" % ARGS.aligner)
  
  BWA.align_standalone(filepaths=ARGS.files,
                       genome  = ARGS.genome,
                       nocleanup = ARGS.nocleanup,
                       nocc = ARGS.nocc)
  


