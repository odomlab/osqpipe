#!/usr/bin/env python
#
# $Id$

'''Script to query a given flowcell ID for lanes of interest, download
the sequencing data (i.e. fastq file) and demultiplex if necessary.'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

# All script code moved into our main pipeline library namespace.
from osqpipe.pipeline.flowcell import FlowCellProcess

###############################################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Process a FlowCell from the LIMS.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Turn on test mode (no filesystem changes).')

  PARSER.add_argument('-n', '--nocheck', dest='checkForLibInDB',
                      action='store_false',
                      help='Deactivate the check for library in repository DB.')

  PARSER.add_argument('flowCell', metavar='<flowCellID>', type=str,
                      help='The flow cell ID (required).')

  PARSER.add_argument('flowLane', metavar='<flowLaneID>', type=int, nargs='?',
                      help='The flow lane ID (optional).')

  PARSER.add_argument('-f', '--force-primary', dest='forcePrimary',
                      action='store_true',
                      help='Start processing a flow cell which is only PRIMARY COMPLETE,'
                      + ' not SECONDARY COMPLETE. It is the user\'s responsibility to ensure'
                      + ' that the LIMS files are ready to be copied across.')
    
  PARSER.add_argument('--force-all', dest='forceAll',
                      action='store_true',
                      help='Start processing a flow cell which is still INCOMPLETE.'
                      + ' Once again, it is the user\'s responsibility to ensure'
                      + ' that the LIMS files are ready to be copied across.')

  PARSER.add_argument('--force-download', dest='forceDownload',
                      action='store_true',
                      help='Force download even if flowcell has already been processed.')

  PARSER.add_argument('-d', '--destination-dir', dest='destdir', type=str,
                      help='Select a different destination directory for downloads'
                      + ' (default is the configured incoming directory).')

  PARSER.add_argument('--trust-lims-adapters', dest='trustAdapt', type=str, default=None,
                      help='If set to the name of a protocol (taken from the Adapter'
                      + ' table) then any empty adapter fields in the database'
                      + ' will be automatically filled using LIMS data. This'
                      + ' bypasses our usual consistency check so should only'
                      + ' be used when certain. It will not overwrite'
                      + ' previously-entered adapter metadata.')

  ARGS = PARSER.parse_args()

  PROC = FlowCellProcess(test_mode        = ARGS.testMode,
                         db_library_check = ARGS.checkForLibInDB,
                         force_primary    = ARGS.forcePrimary,
                         force_all        = ARGS.forceAll,
                         trust_lims_adapters = ARGS.trustAdapt,
                         force_download   = ARGS.forceDownload)

  PROC.run(ARGS.flowCell, ARGS.flowLane, destdir = ARGS.destdir)

