#!/usr/bin/env python
#
# $Id$

'''Script to add a lane-associated file to the repository.'''

import logging
import os.path
from datetime import date

from osqpipe.pipeline.utilities import parse_repository_filename, checksum_file
from osqpipe.models import Filetype, Lane, Lanefile, Facility

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

###########################################################

class RepoFileHandler(object):

  '''Class which is almost certainly overkill given the limited
  functionality left in this script, post-refactor.'''

  def __init__(self):
    LOGGER.setLevel(logging.DEBUG)

  @staticmethod
  def run(fns):
    '''Main entry point for the class.'''

    # Assume that the first file is representative.
    fname = fns[0]
    LOGGER.info(fname)
    (code, facility, lanenum, pipeline) = parse_repository_filename(fname)
    lanelist = Lane.objects.filter(library__code=code,
                                   lanenum=lanenum,
                                   facility__code=facility)
    if len(lanelist) == 0:
      LOGGER.error("Could not find lane for '%s'", fname)
    elif len(lanelist) > 1:
      LOGGER.error("Found multiple lanes for '%s': %s",
                   fname, ", ".join([x.id for x in lanelist]))
    else:
      lane = lanelist[0]
      for fname in fns:
        chksum = checksum_file(fname)
        filetype = Filetype.objects.guess_type(fname)
        basefn = os.path.split(fname)[1]
        LOGGER.debug("basefn: '%s'" % (basefn))
        fnparts = os.path.splitext(basefn)
        if fnparts[1] == '.gz':
          basefn = fnparts[0]
        LOGGER.debug("basefn: '%s'" % (basefn))
        lanefile = Lanefile(filename=basefn, checksum=chksum,
                            filetype=filetype, lane=lane,
                            description='', pipeline=pipeline)
        lanefile.save()

###########################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Add a list of files associated with'
              + ' a single lane to the repository.')

  PARSER.add_argument('files', metavar='<files>', type=str, nargs='+',
                      help='The list of files.')

  ARGS = PARSER.parse_args()

  HND = RepoFileHandler()

  HND.run(ARGS.files)

