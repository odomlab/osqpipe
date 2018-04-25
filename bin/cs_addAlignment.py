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

'''Script to take a list of bed, bam, bgr and wiggle files, record
them in the database and move the files to the correct part of the
repostory filesystem.'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

# The bulk of the code which used to be in this script is now moved
# into our library namespace.
from osqpipe.pipeline.alignment import AlignmentHandler

###########################################################


if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Insert alignment data into the repository. Note that the'
    + ' script currently assumes that all files are gzipped.')

  PARSER.add_argument('files', metavar='<alignment files>', type=str, nargs='+',
                      help='A list of files to be associated with this'
                      + ' alignment. At least one of these must be a bed file.')

  PARSER.add_argument('-p', '--prog', dest='prog', type=str, default='',
                      help='The program associated with these files.')

  PARSER.add_argument('-q', '--params', dest='params', type=str, default='',
                      help='The program parameters.')

  PARSER.add_argument('-g', '--genome', dest='genome', type=str, default='',
                      help='The genome used in the alignment.')

  PARSER.add_argument('-ht', '--headtrim', dest='headtrim', type=int, default=0,
                      help='Whether the tailtrim option was used.')

  PARSER.add_argument('-tt', '--tailtrim', dest='tailtrim', type=int, default=0,
                      help='Whether the headtrim option was used.')

  ARGS = PARSER.parse_args()

  HND = AlignmentHandler(params   = ARGS.params,
                         prog     = ARGS.prog,
                         genome   = ARGS.genome,
                         headtrim = ARGS.headtrim,
                         tailtrim = ARGS.tailtrim)
  
  HND.add(ARGS.files)

