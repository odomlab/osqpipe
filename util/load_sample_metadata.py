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

'''
Simple script to load Characteristic data for Samples in the
repository. The input TSV file should contain headers identifying the
Characteristic name. The first column should contain Sample names. It
is assumed that all Samples and Characteristics are already available
in the repository; this script merely makes links. Pre-existing
Sample-Characteristic links of the designated category will be overwritten
for updates. Any error will result in rollback of the transaction.
'''

import sys
import re

from osqpipe.models import Sample, Characteristic, SizeUnit

import django
django.setup()

from django.db import transaction

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# The file field delimiter (TSV).
DELIMITER = '\t'

@transaction.atomic
def load_metadata_file(infile, sizecat='Diameter', sizeunit='mm', commentcat='Comments'):
  with open(infile, 'r') as infh:
    header = [ field.strip() for field in infh.next().split(DELIMITER) ]
    for line in infh:
      row    = [ field.strip() for field in line.split(DELIMITER) ]
      sample = Sample.objects.get(name__iexact=row[0])
      chars  = dict(zip(header[1:], row[1:]))
      for charname, charval in chars.iteritems():
        if charval == '':
          continue
        LOGGER.debug("Processing %s: %s...", charname, charval)
        if charname.lower() == sizecat.lower():
          sizeunit = SizeUnit.objects.get(name=sizeunit)
          sample.size = charval
          sample.size_unit = sizeunit
          sample.save()
        elif charname.lower() == commentcat.lower():
          if sample.comment is None:
            # new comment
            note = charval
          elif re.search(re.escape(charval), sample.comment):
            # comment already present
            continue
          else:
            # append comment
            note = " ".join((sample.comment, charval))
          sample.comment = note
          sample.save()
        else:

          # Handle pre-existing Characteristics in the charval
          # category. This is potentially destructive if we ever move
          # to a data model where a Sample can have more than one
          # Characteristic of a given category.
          found = False
          preexisting = Characteristic.objects.filter(category__iexact=charname,
                                                      samples=sample)
          for char in preexisting:
            if char.value.lower() == charval.lower():
              LOGGER.debug("Sample %s already linked to Characteristic [%s]", sample, char)
              found = True
            else:
              LOGGER.warning("Removing mismatched Characteristic for Sample %s: [%s]", sample, char)
              sample.characteristics.remove(char)
          if found:
            continue

          # If an appropriate Characteristic is not already linked,
          # create a new link.
          try:
            char = Characteristic.objects.get(category__iexact=charname,
                                              value__iexact=charval)
          except Characteristic.DoesNotExist:
            LOGGER.error("Unable to find Characteristic in database: [%s: %s]", charname, charval)
            sys.exit(1) # Will rollback transaction.

          LOGGER.info("Creating link between Sample %s and Characteristic [%s]", sample, char)
          sample.characteristics.add(char)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='''Link repository Samples with Characteristics. The input TSV file
should contain headers identifying the Characteristic name. The first
column should contain Sample names. It is assumed that all Samples and
Characteristics are already available in the repository; this script
merely makes links. Any error will result in rollback of the
transaction.''')

  PARSER.add_argument('-f', '--file', dest='infile', type=str, required=True,
                      help='The name of the TSV file containing sample metadata. ')

  PARSER.add_argument('--size-category', dest='sizecat', type=str, default='Diameter',
                      help='The column containing size metadata for the samples.')

  PARSER.add_argument('--size-unit', dest='sizeunit', type=str, default='mm',
                      help='The unit to be used for sample size metadata.')

  PARSER.add_argument('--comment-category', dest='commentcat', type=str, default='Comments',
                      help='The column containing a comment for the samples.')

  ARGS = PARSER.parse_args()

  load_metadata_file(infile     = ARGS.infile,
                     sizecat    = ARGS.sizecat,
                     sizeunit   = ARGS.sizeunit,
                     commentcat = ARGS.commentcat)
