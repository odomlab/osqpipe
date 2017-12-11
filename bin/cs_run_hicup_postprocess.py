#!/usr/bin/env python
#
# $Id$
#
import os
import sys

if __name__ == '__main__':
    
    import argparse

    PARSER = argparse.ArgumentParser(
        description='Runs HiCUP results post-processing in the cluster. Use together with cs_run_hicup.py.')
    PARSER.add_argument('--fq1', dest='fq1', type=str,
                        help='Fastq of pair1', required=True)

    ARGS = PARSER.parse_args()


    HC = HiCUP(fq1=ARGS.fq1)
    HC.write_hicup_config()
    HC.run_hicup()
