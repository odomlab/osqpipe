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
Script used to update HCC project bam files such that the sample tag
reflects library.sample.name. Currently in testing; we may wish to
extend the tag to include tumour number as well as mouse ID.
'''

import os
from logging import INFO
from osqutil.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Alnfile, Library
from osqutil.utilities import checksum_file, call_subprocess, sanitize_samplename
from osqutil.config import Config
from osqpipe.pipeline.gatk import update_bam_readgroups

from django.db import transaction
from shutil import move

CONFIG = Config()

@transaction.atomic
def update_repo_bamfile(bam, newfile):
  '''
  Carefully replace a bam file in the repository.
  '''
  # Ensure we're using the latest bam object from the database within
  # this transaction.
  bam = Alnfile.objects.get(id=bam.id)
  checksum = checksum_file(newfile, unzip=False)
  deleteme = "%s.bak" % (bam.repository_file_path,)
  move(bam.repository_file_path, deleteme)
  move(newfile, bam.repository_file_path)
  bam.checksum = checksum
  bam.save()
  os.unlink(deleteme)

def update_library_bam_readgroups(libcode, update=False):
  '''
  Add or Replace the read groups for all the bam files attached to a
  given library. If update is False (the default), picard
  AddOrReplaceReadGroups is used; if update is True then just the bam
  file header is rewritten using an internal pipeline function.
  '''
  lib  = Library.objects.get(code=libcode)
  bams = Alnfile.objects.filter(alignment__lane__library__code=libcode,
                                filetype__code='bam')
  common_args = ('VALIDATION_STRINGENCY=SILENT',
                 'TMP_DIR=%s' % CONFIG.tmpdir)

  for bam in bams:
    LOGGER.info("Updating bam file: %s", bam.filename)
    checksum = checksum_file(bam.repository_file_path, unzip=False)
    if checksum != bam.checksum:
      raise ValueError("Stored bam checksum does not agree with that in the repository.")
    if update:
      LOGGER.debug("Rewriting bam file header: %s", bam.filename)
      update_bam_readgroups(bam)
    else:
      tmpfile = "%s.update_rg" % (bam.filename,)
      cmd = ('picard', 'AddOrReplaceReadGroups',
             'INPUT=%s'  % bam.repository_file_path,
             'OUTPUT=%s' % tmpfile,
             'RGLB=%s'   % lib.code,
             'RGSM=%s'   % sanitize_samplename(lib.sample.name),
             'RGCN=%s'   % bam.alignment.lane.facility.code,
             'RGPU=%d'   % int(bam.alignment.lane.lanenum),
             'RGPL=illumina') + common_args

      LOGGER.debug("Running command: %s", " ".join(cmd))
      call_subprocess(cmd, path=os.environ['PATH'])
      update_repo_bamfile(bam, tmpfile)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description=\
       'Update bam file read groups and repository checksums to match'
                          + ' current HCC project best practice.')

  PARSER.add_argument('-l', '--library', dest='library', type=str, required=True,
                      help='The library code to process.')

  PARSER.add_argument('-u', '--update', dest='update', action='store_true',
                      help='Update mode only (requires the bam files to already'
                      + ' have had read groups added). This is faster than the'
                      + ' default alternative.')

  ARGS = PARSER.parse_args()

  update_library_bam_readgroups(ARGS.library, update=ARGS.update)
