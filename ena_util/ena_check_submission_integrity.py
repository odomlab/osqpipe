#!/usr/bin/env python
#
# $id$

import sys
import os

# set up logger
from osqutil.setup_logs import configure_logging
from logging import WARNING
LOGGER = configure_logging(level=WARNING)

# import config
from osqutil.config import Config

# For insertion of lane info:
import django
from osqpipe.models import Lane, Library, ExternalRecord

# set up config
DBCONF = Config()

django.setup()

def check_ena_submission_integrity(code):

    library = None
    try:
        library = Library.objects.get(code=code)
    except Library.DoesNotExist:
        LOGGER.error("Library with code=%s not found!", code)
        sys.exit(1)

    # check for external ENA record for library
    try:
        ers = ExternalRecord.objects.filter(libraries=library, repository__name='ENA')
    except ExternalRecord.DoesNotExist:
        LOGGER.error("No external records associated with library %s.", code)
        sys.exit(1)

    # check for external ENA record for library associated sample
    try:
        ers = ExternalRecord.objects.filter(samples=library.sample, repository__name='ENA')
    except ExternalRecord.DoesNotExist:
        LOGGER.error("No external records associated with sample %s.", self.obj)
        sys.exit(1)

    # check for external ENA record for each of the lanes
    try:
        lanes = Lane.objects.filter(library__code=code)
    except Lane.DoesNotExist:
        LOGGER.error("No lanes for library %s.", code)
        sys.exit(1)
    if len(lanes) == 0:
        LOGGER.error("No lanes for library %s.", code)
        sys.exit(1)
    for lane in lanes:
        try:
            ers = ExternalRecord.objects.filter(lanes=lane, repository__name='ENA')
        except ExternalRecord.DoesNotExist:
            LOGGER.error("No external records associated with lane %s.", lane.id)
            sys.exit(1)
    print "%s\tSubmitted" % code

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
      description='Check ENA submission integrity for a library. The script checks for association of ENA for library and related sample and lanes.')

  PARSER.add_argument('-l', '--library', dest='library', type=str,
                      help='Library code (donumber). E.g. do1234.')

  ARGS = PARSER.parse_args()

  check_ena_submission_integrity(ARGS.library)
