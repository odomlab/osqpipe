#!/usr/bin/env python
#
# $Id$

'''Script to scan through all alignments in the repository and create
bedgraph files for those which need them.'''

import os
import os.path
import logging
from datetime import date
from pipes import quote

from django.db import transaction

from osqpipe.pipeline.utilities import call_subprocess, checksum_file
from osqpipe.models import Filetype, Library, Lane, Alignment, Alnfile, Facility
from osqpipe.pipeline.config import Config

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

BED2BGR = "makeWiggle -B -1 %s %s"

###############################################################################

class BedGraphCreator(object):

  '''Core class used to iterate over all alignments in the database.'''

  __slots__ = ('testMode', 'conf', 'bedtype', 'bgrtype')

  def __init__(self, testMode=False):
    LOGGER.setLevel(logging.INFO)
    self.testMode = testMode
    self.conf     = Config()
    self.bedtype  = Filetype.objects.get(code='bed')
    self.bgrtype  = Filetype.objects.get(code='bgr')

  def needs_bed_graph(self, aln):
    '''Look at the files associated with the alignment; if there is a
    bed file but no bgr file then return true, else return false.'''
    bgr = aln.alnfile_set.filter(filetype=self.bgrtype)
    bed = aln.alnfile_set.filter(filetype=self.bedtype).exclude(filename__contains='chr21')
    return (len(bed) == 1 and len(bgr) == 0)

  def load_library(self, aln):
    '''Convenience function to load the lane and library for a given
    alignment.'''
    lane = aln.lane
    return (lane, lane.library, lane.facility.code)

  @transaction.commit_on_success
  def make_bed_graph(self, aln):
    '''Code wrapper for makeWiggle.'''
    bed = aln.alnfile_set.filter(filetype=self.bedtype).exclude(filename__contains='chr21')[0]
    
    # Note that this doesn't use alnfile.repository_file_path, because
    # it's allegedly working directly with uncompressed BED and BGR
    # files stored in the repository. Since it's not respecting the
    # filetype.gzip flags, I surmise that this script is old and not
    # being used at the moment. I've therefore not gone to the trouble
    # of writing scads of gzip-handling code to fix it at this point
    # (TFR).
    lib   = aln.lane.library
    bedFN = os.path.join(self.conf.repositorydir,
                         lib.code,
                         bed.filename)
    bgrBASE = os.path.join(self.conf.repositorydir,
                           lib.code,
                           os.path.splitext(bed.filename)[0])
    bgrFN = bgrBASE + self.bgrtype.suffix
    cmd = BED2BGR % (quote(bedFN), quote(bgrBASE))
    LOGGER.debug(cmd)
    if not self.testMode:
      call_subprocess(cmd, shell=True, path=self.conf.hostpath)
      if not os.path.exists(bgrFN):
        LOGGER.error("Failed to create bgr file '%s'" % (bgrFN,))
      else:
        chksum = checksum_file(bgrFN)
        bgr = Alnfile(filename=os.path.basename(bgrFN), checksum=chksum,
                      filetype_id=self.bgrtype.id, description='',
                      alignment=aln)
        bgr.save()

  def run(self):

    '''Main entry point for the class.'''

    for aln in Alignment.objects.all():
      if self.needs_bed_graph(aln):
        (lane, lib, facility) = self.load_library(aln)
        LOGGER.info("%s %s%02d", lib.code, facility, lane.lanenum)
        self.make_bed_graph(aln)
      else:
        (lane, lib, facility) = self.load_library(aln)
        LOGGER.debug("skipping %s %s%02d",
                     lib.code, facility, lane.lanenum)

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

