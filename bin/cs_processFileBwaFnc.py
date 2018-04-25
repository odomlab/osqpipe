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

'''Script to align fastq files against their genome as specified in
the repository, using bwa.'''

import sys

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.pipeline.file_processor import FileProcessingManager

###############################################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Submit files for processing on the cluster.')

  PARSER.add_argument('files', metavar='<filenames>', type=str, nargs='+',
                      help='The name of the file or files to process.\n'
                      + 'The files may be a single or a pair of fastq files; or, when used together with option --facility=\'SAN\', a bam file.')

  PARSER.add_argument('-o', '--opt', dest='options', type=str,
                      action='append', default=[],
                      help='Options to be passed through;'
                      + ' e.g., trimtail=36 (trims sequence to 36b), noalign (skip aligning), convert (force conversion of quality scores from Solexa 1.5 to Sanger format).')

  PARSER.add_argument('-f', '--facility', dest='facility',
                      type=str, default='CRI',
                      help='The facility supplying the input'
                      + ' file (e.g. SAN for Sanger).'
                      + ' In case input is a bam file and --facility=\'SAN\','
                      + ' the script expects presence of .bam.meta file containing metadata from Sanger IRODs.\n')

  PARSER.add_argument('--aligner', type=str, dest='aligner', choices=('bwa', 'tophat', 'star'), default=None, 
                      help='Tophat is the default aligner for rnaseq data while bwa is default for all other library types.')
  
  PARSER.add_argument('--algorithm', type=str, dest='algorithm', choices=('aln', 'mem'),
                      help='The bwa algorithm to use (aln or mem). The default behaviour'
                      + ' is to pick the algorithm based on the read length in the fastq files.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Turn on test mode.')

  PARSER.add_argument('--paired', dest='paired', action='store_true',
                      help='Force bam2fastq export to treat the data files as paired end (Sanger pipeline).')

  ARGS = PARSER.parse_args()

  if not (len(ARGS.files) in (1, 2)):
    sys.exit('Either one or two filenames must be passed to this script.')

  # Parse out our --opt arguments into a dict
  OPTIONS = {}
  for a in ARGS.options:
    opt_fields = a.split("=")
    opt_fields = [ opt.strip() for opt in opt_fields ]
    if len(opt_fields) == 1:
      OPTIONS[opt_fields[0]] = True
    else:
      val = "=".join(opt_fields[1:])
      OPTIONS[opt_fields[0]] = val

  # Best to use the TEST_MODE global here.
  TEST_MODE = ARGS.testMode

  FPM = FileProcessingManager(options=OPTIONS,
                              facility=ARGS.facility,
                              force_paired_end=ARGS.paired,
                              bwa_algorithm=ARGS.algorithm,
                              test_mode=ARGS.testMode, aligner=ARGS.aligner)
  FPM.run(ARGS.files)
