#!/usr/bin/env python

'''
Script which bridges the gap between our standard sequencing pipeline
and the Bioinformatics Core GATK preprocessing pipeline.

The steps included in this script are:

1. Merge all the files from the specified libraries (which should all
come from the same HCC nodule).

2. Transfer the merged bam file to the cluster.

2. Run picard MarkDuplicates on the merged bam file.

3. Run picard BuildBamIndex on the output.

4. Start the GATK IndelRealigner-BaseRecalibrator pipeline as provided
by the Bioinformatics Core pipeline.

5. Submit cleanup jobs to transfer the output back to local host, and
delete working files on the cluster. The GATK log files are currently
retained in the working directory.
'''

from osqpipe.pipeline.gatk import GATKPreprocessor

from logging import INFO
from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description=\
               'Script to initiate the HCC GATK preprocessing pipeline.')
  
  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=False,
                      help='The alignment genome used to filter the input files.')

  PARSER.add_argument('-l', '--libtype', dest='libtype', type=str, required=False,
                      help='The library type (genome, exome) used to filter the input files.')

  GROUP = PARSER.add_mutually_exclusive_group(required=False)

  # I'd personally love to have this on GROUP, but it's not compatible
  # with nargs='*' and insists on nargs='?'. I personally think this
  # is a bug in argparse.
  PARSER.add_argument('libraries', metavar='<libcodes>', type=str, nargs='*',
                      help='The names of the libraries to merge and load'
                      + ' into the pipeline. All the files on the command'
                      + ' line should come from the same HCC nodule. Either this,'
                      + ' the --sample, '
                      + ' or the --merged-bam arguments must be supplied.')

  GROUP.add_argument('-s', '--sample', dest='sample', type=str,
                     help='The Sample ID to process; all libraries linked'
                     + ' to this sample (optionally filtered by genome and library type)'
                     + ' will be combined.')

  GROUP.add_argument('-m', '--merged-bam', dest='mergedbam', type=str,
                     help='The name of a merged bam file to use instead'
                     + ' of library codes or sample ID.')

  ARGS = PARSER.parse_args()

  PROC = GATKPreprocessor()

  if ARGS.mergedbam is not None:
    PROC.gatk_preprocess_bam(ARGS.mergedbam)
  elif ARGS.sample is not None:
    PROC.gatk_preprocess_sample(ARGS.sample, genome=ARGS.genome, libtype=ARGS.libtype)
  else:
    PROC.gatk_preprocess_libraries(ARGS.libraries, genome=ARGS.genome)
