#!/usr/bin/env python
#
# $Id$
#
import os
import sys

from osqpipe.pipeline.hicup import HiCUP

if __name__ == '__main__':
    
    import argparse

    PARSER = argparse.ArgumentParser(
        description='Runs HiCUP on read file or file pair. Note that the program leaves results in cluster. Use together with cs_run_hicup_postprocess.py.')
    PARSER.add_argument('--fq1', dest='fq1', type=str,
                        help='Fastq of pair1', required=True)
    PARSER.add_argument('--fq2', dest='fq2', type=str, default=None,
                        help='Fastq of pair2')
    PARSER.add_argument('--genome', dest='genome', type=str, required=True,
                     help='Reference genome name.')
    PARSER.add_argument('--enzyme', dest='enzyme', type=str, required=True,
                        help='Name of the restriction enzyme which was used in the HiC assay.')

    ARGS = PARSER.parse_args()


    HC = HiCUP(fq1=ARGS.fq1, genome=ARGS.genome, enzyme=ARGS.enzyme, fq2=ARGS.fq2)
    HC.write_hicup_config()
    HC.run_hicup()
