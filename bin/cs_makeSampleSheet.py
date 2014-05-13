#!/usr/bin/env python
#
# $Id$
#
'''This is currently a very underdeveloped script that's been used to
partially generate the sample sheets used with the demuxIllumina
program. Note that the output will almost invariably need editing; the
main benefit of running this script is its automated handling of
adapter sequences.'''

from logging import INFO

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

from osqpipe.pipeline.utilities import build_incoming_fastq_name

from osqpipe.models import Library

class DemuxSheetMaker(object):

  '''Very simple wrapper class used to select adapter sequences and
  output them to a sample sheet. Note that this is vastly
  overengineered for its current function.'''

  __slots__ = ('verbose')

  def __init__(self, verbose=False):
    self.verbose = verbose
    LOGGER.setLevel(INFO)

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
      fdesc.write("%s %s\n" % (lib.adapter.sequence, fname))

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
