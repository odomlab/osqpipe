#!/usr/bin/env python
#
# $Id$

""" Merge BWA generated files """
__author__ = "Margus Lukk"
__date__ = "05 Mar 2012"
__version__ = "0.2"
__credits__ = "Adjusted from cs_runMaqWithSplit_Merge written by Gordon Brown."

# Known bugs: 
# 

import sys
import os
import os.path
import re

from osqpipe.pipeline.setup_logs import configure_logging
from logging import INFO, WARNING
LOGGER = configure_logging(level=INFO)

from osqpipe.pipeline.bwa_runner import BwaAlignmentManager

################################################################################

      
###############   M A I N   P R O G R A M   ############

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Merge a set of aligned BAM files into one.')

  PARSER.add_argument('outfile', metavar='<output BAM filename>', type=str,
                      help='The name of the output BAM file.')

  PARSER.add_argument('infiles', metavar='<input BAM file(s)>',
                      type=str, nargs='+',
                      help='The BAM files to merge.')

  PARSER.add_argument('--loglevel', type=int, dest='loglevel', default=WARNING,
                      help='The level of logging.')

  PARSER.add_argument('--rcp', type=str, dest='rcp',
                      help='Remote file copy (rcp) target.')

#  Currently unused:  
#  PARSER.add_argument('--force', dest='force', action='store_true',
#                      default=False,
#                      help='Overwrite existing output file.')

  PARSER.add_argument('--cleanup', dest='cleanup', action='store_true',
                      help='Delete all temporary files.')

  PARSER.add_argument('--group', type=str, dest='group',
                      help='The user group for the files.')

  PARSER.add_argument('-d', '--debug', dest='debug', action='store_true',
                      help='Turn on debugging output.')

  ARGS = PARSER.parse_args()

  BSUB = BwaAlignmentManager(debug      = ARGS.debug,
                             cleanup    = ARGS.cleanup,
                             loglevel   = ARGS.loglevel,
                             group      = ARGS.group)

  BSUB.merge_alignments(input_fns  = ARGS.infiles,
                        output_fn  = ARGS.outfile,
                        rcp_target = ARGS.rcp)

  
