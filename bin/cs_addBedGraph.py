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
Script to scan through all alignments in the repository and create
bedgraph files for those which need them.
'''

import os
import os.path
from datetime import date
from pipes import quote
from shutil import move

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from django.db import transaction

from osqutil.utilities import call_subprocess, checksum_file, rezip_file, set_file_permissions
from osqpipe.models import Filetype, Library, Lane, Alignment, Alnfile, Facility
from osqutil.config import Config

BED2BGR = "makeWiggle -B -1 %s %s"

###############################################################################

class BedGraphCreator(object):
  '''
  Core class used to iterate over all alignments in the database.
  '''
  __slots__ = ('testMode', 'conf', 'bedtype', 'bgrtype')

  def __init__(self, testMode=False):
    self.testMode = testMode
    self.conf     = Config()
    self.bedtype  = Filetype.objects.get(code='bed')
    self.bgrtype  = Filetype.objects.get(code='bgr')

  def needs_bed_graph(self, aln):
    '''
    Look at the files associated with the alignment; if there is a
    bed file but no bgr file then return true, else return false.
    '''
    bgr = aln.alnfile_set.filter(filetype=self.bgrtype)
    bed = aln.alnfile_set.filter(filetype=self.bedtype).exclude(filename__contains='chr21')
    return (len(bed) == 1 and len(bgr) == 0)

  @transaction.atomic
  def make_bed_graph(self, aln):
    '''
    Code wrapper for makeWiggle.
    '''
    aln = Alignment.objects.get(id=aln.id) # Reload passed object within transaction.
    bed = aln.alnfile_set.filter(filetype=self.bedtype).exclude(filename__contains='chr21')[0]
    
    # Note makeWiggle can read gzipped bed files directly; we use that fact here.
    lib   = aln.lane.library
    bedFN = bed.repository_file_path

    # Write to local directory first.
    bgrBASE = os.path.splitext(bed.filename)[0]
    bgrFN   = bgrBASE + self.bgrtype.suffix
    cmd = BED2BGR % (quote(bedFN), quote(bgrBASE))
    LOGGER.debug(cmd)
    if not self.testMode:
      call_subprocess(cmd, shell=True, path=self.conf.hostpath)
      if not os.path.exists(bgrFN):
        LOGGER.error("Failed to create bgr file '%s'" % (bgrFN,))
      else:
        chksum = checksum_file(bgrFN)
        bgr = Alnfile(filename=os.path.basename(bgrFN), checksum=chksum,
                      filetype=self.bgrtype,
                      alignment=aln)
        bgrFN = rezip_file(bgrFN)
        move(bgrFN, bgr.repository_file_path)
        set_file_permissions(self.conf.group, bgr.repository_file_path)
        bgr.save()

  def run(self):
    '''
    Main entry point for the class.
    '''
    for aln in Alignment.objects.all():
      if self.needs_bed_graph(aln):
        LOGGER.info("Creating bedgraph file for alignment %s", aln)
        self.make_bed_graph(aln)
      else:
        LOGGER.debug("Skipping alignment %s", aln)

###############################################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Add BedGraph data to the repository'
              + ' for any alignment which needs them.')

  PARSER.add_argument('-t', '--test', dest='testMode', action='store_true',
                      help='Turn on test mode.')

  ARGS = PARSER.parse_args()

  CREATOR = BedGraphCreator(testMode=ARGS.testMode)

  CREATOR.run()

