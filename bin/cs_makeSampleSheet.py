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

'''This is currently a very underdeveloped script that's been used to
partially generate the sample sheets used with the demuxIllumina
program. Note that the output will almost invariably need editing; the
main benefit of running this script is its automated handling of
adapter sequences.'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqutil.utilities import build_incoming_fastq_name

from osqpipe.models import Library

class DemuxSheetMaker(object):

  '''Very simple wrapper class used to select adapter sequences and
  output them to a sample sheet. Note that this is vastly
  overengineered for its current function.'''

  __slots__ = ('verbose')

  def __init__(self, verbose=False):
    self.verbose = verbose

  def run(self, flowCellID, lane, libs, outfile):
    '''The main entry point for the class.'''
        
    fdesc = open(outfile, 'w')
    for code in libs:
      try:
        lib = Library.objects.get(code=code)
      except Library.DoesNotExist, err:
        raise StandardError("Unable to find library code %s in the database." % code)
      if not lib.adapter:
        raise StandardError("Library has no adapter associated: %s." % code)
      fname = build_incoming_fastq_name(code.lower(), flowCellID, lane, 1)
      if lib.adapter2 is None:
        fdesc.write("%s\t%s\t%s\n" % (lib.adapter.sequence, fname, lib.code))
      else:
        fdesc.write("%s\t%s\t%s\t%s\n" % (lib.adapter.sequence, lib.adapter2.sequence, fname, lib.code))

    fdesc.close()

if __name__ == '__main__':

  import argparse
  
  PARSER = argparse.ArgumentParser(
    description='Somewhat kludgey script to create a sample sheet'
    + ' for use with demuxIllumina.')

  PARSER.add_argument('flowCell', metavar='<flowCellID>', type=str,
                      help='The flow cell ID (required).')

  PARSER.add_argument('flowLane', metavar='<flowLaneID>', type=str,
                      help='The flow lane ID (required). Note that'
                      + ' for handling data from multiple lanes you'
                      + ' will have to manually edit the output sheet.')

  PARSER.add_argument('libraries', metavar='<library IDs>', type=str, nargs='+',
                      help='A list of library IDs.')

  PARSER.add_argument('-o', '--out', dest='outfile', type=str, required=True,
                      help='The name of the output file.')

  PARSER.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                      help='Turn on verbose output.')

  ARGS = PARSER.parse_args()

  DSM = DemuxSheetMaker(verbose=ARGS.verbose)

  DSM.run(flowCellID = ARGS.flowCell,
          lane       = ARGS.flowLane,
          libs       = ARGS.libraries,
          outfile    = ARGS.outfile)
