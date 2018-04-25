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

'''Script to fetch fastq files from the LIMS.'''

import sys
import os
import os.path
import re

from osqutil.utilities import build_incoming_fastq_name, unzip_file, \
    set_file_permissions, checksum_file
from .upstream_lims import Lims
from ..models import Library, Lane
from osqutil.config import Config

from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging('fetch_fastq')

###############################################################################

class FQFileFetcher(object):

  '''Class used to query the LIMS for fastq files associated with a
  given flowcell and download them to a destination directory.'''

  __slots__ = ('destination', 'lims', 'targets', 'test_mode', 'conf', 'unprocessed_only','force_download')

  def __init__(self, destination, lims=None, test_mode=False, unprocessed_only=False, force_download=False):

    self.conf        = Config()
    self.test_mode   = test_mode
    self.unprocessed_only = unprocessed_only
    self.destination = destination
    self.force_download = force_download
    self.targets = set()    
    if lims is None:
      lims = Lims()
    if not lims.running():
      LOGGER.error("Remote LIMS access broken... cannot continue.")
      sys.exit("LIMS not running.")
    self.lims    = lims

    if self.test_mode:
      LOGGER.setLevel(DEBUG)
    else:
      LOGGER.setLevel(INFO)

  def retrieve_fqfile(self, lfile, libname):
    '''
    Given a LimsFile object and a library name, retrieve the actual
    fastq files from the LIMS and store it in self.destination.
    '''
    if not os.path.exists(self.destination):
      LOGGER.error("Destination '%s' does not exist.", self.destination)
      return

    filename = lfile.uri.split('/')[-1]
    LOGGER.info("LIMS File: %s", filename)

    # Current file naming convention:
    # <SLX ID>.<run number>.s_<lane no>.r_<N>.fq.gz
    # Where N=1 indicates single end sequencing
    #       N=2 indicates paired end
    #       N=3 indicates paired end, multiplexed
    #       N=4 indicates dual-indexed.
    #
    # Old file naming convention (still supported):
    # s_<lane no>(?:_<N>)_sequence.txt.gz
    # Where N is only present for non-single-end sequencing
    # and means the same as above.

    # This regex supports both old and new naming conventions.
    fnpat = re.compile( # FIXME see parse_incoming_fastq_name
      r"SLX\-\d+\.[\.\w-]+\.s_\d+\.r_(\d+).fq.gz$"
      + r"|s_\d+(?:_(\d+))?_sequence.txt.gz$"
      + r"|s_\d+\.(tar)$")  # 10X_FASTQ_TAR
    matchobj = fnpat.search(filename)

    if matchobj is None:
      LOGGER.error("FASTQ file does not conform to"
                   + " known naming convention: %s", filename)
      return

    flowpair = 1
    fastqtar = False
    if matchobj.group(1) is not None:  # New naming
      stype = int(matchobj.group(1))
      if stype > 1: # it's a paired-end-style name
        flowpair = stype
    elif matchobj.group(2) is not None: # Old naming
      flowpair = int(matchobj.group(1))
    elif matchobj.group(3) == 'tar':   # 10X FASTQ tar file, probably.
      fastqtar = True
    else:
      LOGGER.error("FASTQ file name regex gave unexpected"
                   + " results; probable error in regex code.")
      return

    # following if statement has been added by Margus to deal
    # properly PE multiplexed lanes where second read suffixed _3.fq
    # rather than _2.fq.
    if (flowpair == 3):
      flowpair = 2
    elif (flowpair == 4):
      LOGGER.warning("Dual indexed files not yet supported. However, in current"
                     + " usage this is likely to be mislabeled by LIMS as part"
                     + " of a wider flowcell annotation.")
#        return   # If we ever start dual indexing this will need to change
      flowpair = 2

    sample_id = libname.lower().replace(" ", "") # do not use underscores here.
    dst = build_incoming_fastq_name(sample_id,
                                    lfile.lane.flowcell.fcid,
                                    lfile.lane.lane,
                                    flowpair)

    if fastqtar: # s/.fq$/.tar/
      dst = os.path.splitext(dst)[0] + '.tar'

    # If the final file has already been downloaded we skip the download.
    target = os.path.join(self.destination, dst)
    if os.path.exists(target):
      LOGGER.warning("Destination file '%s' exists. Cannot overwrite.", target)
      return

    compressed = False
    if os.path.splitext(lfile.uri)[1] == self.conf.gzsuffix:
      compressed = True
      dst += self.conf.gzsuffix
    target = os.path.join(self.destination, dst)

    # We also refuse to download over an intermediary gzipped
    # file. This is more likely to be an error so we raise an
    # exception here. If the file is good it should have been
    # uncompressed already.
    if os.path.exists(target):
      raise StandardError("Download location '%s' exists. Cannot overwrite." % target)

    # We download these over http now.
    LOGGER.debug("Downloading LIMS file ID %s to %s", lfile.lims_id, target)
    if not self.test_mode:

      # This is actually the preferred download mechanism.
      try:
        self.lims.get_file_by_uri(lfile.uri, target)

      # Fall back to download via LIMS API, if supported.
      except Exception, err:
        if lfile.lims_id is not None:
          self.lims.get_file_by_id(lfile.lims_id, target)
        else:
          raise err

      set_file_permissions(self.conf.group, target)

    if not os.path.exists(target) and not self.test_mode:
      LOGGER.error("Failed to retrieve file '%s'", dst)
    else:

      # Compare the md5sum against those available in upstream LIMS.
      if not self.test_mode and lfile.md5sum is not None:
        md5 = checksum_file(target, unzip=False)
        if md5 != lfile.md5sum:
          raise StandardError("File md5sum (%s) disagrees with upstream LIMS (%s): %s"
                              % (md5, lfile.md5sum, target))

    # Files are typically still compressed at this stage. This should
    # be handled seamlessly by downstream code.
    self.targets.add(target)
    LOGGER.info("Downloaded file to %s", target)


  def retrieve_fqfiles(self, lane, libname):
    '''Given a Lane object and a library name, retrieve fastq files
    from the LIMS itself. This method will preferentially retrieve
    demultiplexed sample FASTQ files first, and only retrieve
    lane-based FASTQ if demultiplexed data are not yet available.'''

    if not os.path.exists(self.destination):
      LOGGER.error("Destination '%s' does not exist.", self.destination)
      return
    files = []
    samples = lane.lims_demuxed_samples()
    if len(samples) > 0: # Lane has demultiplexed files available for
                         # at least one sample.

      # We check here that all sample demultiplexed files are
      # available.
      expected = lane.lims_samples()
      if not all([ x in samples for x in expected ]):
        raise StandardError("LIMS has demultiplexed files for only some"
                            + " of the expected lane samples.")

      # Deal with demultiplexed files
      LOGGER.info("Downloading demultiplexed fastq files per sample.")
      for libname in samples:

        if Lane.objects.filter(flowcell=lane.flowcell.fcid, flowlane=lane.lane, library__code__iexact=libname).exists():
          if self.unprocessed_only:
            if not self.force_download:
              LOGGER.info("Skipping files for lane already existing in repository.")
              continue
            else:
              LOGGER.info("Forcing download. Lane already existing in repository.")
          else:
            LOGGER.warning("Downloading files for lane already existing in repository.")

        files = lane.lims_demuxed_files(libname)
        if len(files) == 0:
          LOGGER.info("No Demultiplexed FASTQ files to retrieve for %s_%s",
                      lane.flowcell.fcid, libname)
        for lfile in files:
          self.retrieve_fqfile(lfile, libname)

    else:

      # Download lane-based files.
      LOGGER.info("Downloading fastq files per lane (demultiplexing"
                  + " locally if necessary).")

      if Lane.objects.filter(flowcell=lane.flowcell.fcid, flowlane=lane.lane, library__code__iexact=libname).exists():
        if self.unprocessed_only:
          LOGGER.info("Skipping files for lane already existing in repository.")
          return
        else:
          LOGGER.warning("Downloading files for lane already existing in repository.")

      files = lane.lims_files('FASTQ')
      if len(files) == 0:
        LOGGER.info("No Lane FASTQ files to retrieve for %s_%d",
                    lane.flowcell.fcid, lane.lane)
      for lfile in files:
        self.retrieve_fqfile(lfile, libname)

    return

  def fetch(self, flowcell, flowlane):
    '''Main entry point for the class. Takes a flowcell ID and a lane
    number, and downloads the fastq file to the destination directory
    configured for this FQFileFetcher object.'''

    flowlane = int(flowlane)
    lims_fc = self.lims.load_flowcell(flowcell)
    if self.test_mode:
      lims_fc.dump()
    lims_lane = lims_fc.get_lane(flowlane)
    libname = lims_lane.user_sample_id
    try:
      lib = Library.objects.search_by_name(libname)
      libname = lib.code
    except Library.DoesNotExist, _err:
      LOGGER.debug("Unable to find library with name: %s", libname)

    self.retrieve_fqfiles(lims_lane, libname)

