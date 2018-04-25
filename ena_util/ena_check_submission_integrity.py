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
