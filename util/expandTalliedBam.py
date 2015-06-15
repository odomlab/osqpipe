#!/usr/bin/env python

'''
A script which can be used to expand an Odom Lab pipeline generated
bam file, in which identical small RNA reads are compressed into a
single record, into a more conventional bam file format (one read per
record).
'''

from osqpipe.pipeline.bampy import open_bamfile

################################################################################

if __name__ == '__main__':

  from argparse import ArgumentParser

  P = ArgumentParser(description='Convert reaper/tally output bam files into'
                     + ' conventional bam format (one read per record).')

  P.add_argument('-t', '--tallied', dest='input', type=str, required=True,
                 help='The reaper/tally-generated bam file to expand.')

  P.add_argument('-e', '--expanded', dest='output', type=str, required=True,
                 help='The name of the output expanded bam file.')

  P.add_argument('-m', '--multi', dest='multi', action='store_true',
                 help='Flag indicating that the script should keep any'
                 + ' multi-mapping reads in the final output.')

  ARGS = P.parse_args()

  with open_bamfile(ARGS.input) as bam:
    bam.expand_tallied_reads(ARGS.output, keep_multi=ARGS.multi)
    
