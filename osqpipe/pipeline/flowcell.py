#!/usr/bin/env python
#
# $Id$

'''Script to query a given flowcell ID for lanes of interest, download
the sequencing data (i.e. fastq file) and demultiplex if necessary.'''

import sys
import os
import os.path
import re

from utilities import parse_incoming_fastq_name, checksum_file, \
    build_incoming_fastq_name, call_subprocess, set_file_permissions, \
    munge_cruk_emails
from upstream_lims import Lims
from config import Config
from ..models import Library, Lane, Status, LibraryNameMap, User

from fetch_fastq import FQFileFetcher

from setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging('flowcell')

###############################################################################

def demux_code(code):
  '''Split a comma-separated list of library codes into a list.'''
  return [ x.strip() for x in code.split(",") ]

###############################################################################

class FlowCellProcess(object):

  '''Main class used for processing flowcells.'''

  __slots__ = ('test_mode', 'db_library_check', 'demux_prog', 'conf', 'ready',
               'lims', '_demux_files', 'output_files')

  def __init__(self, test_mode=False,
               db_library_check=True, demux_prog='demuxIllumina',
               force_primary=False, lims=None):

    self.conf             = Config()
    self.test_mode        = test_mode
    self.db_library_check = db_library_check
    self.demux_prog       = demux_prog
    self.ready            = 'COMPLETE'

    # This may now be obsolete with the transition to Genologics LIMS.
    if force_primary:
      self.ready = (self.ready, 'PRIMARY COMPLETE')

    self._demux_files    = {}
    self.output_files    = []
    if lims is None:
      lims = Lims()
    if not lims.running():
      LOGGER.error("Remote LIMS access broken... cannot continue.")
      sys.exit("LIMS not running.")
    self.lims = lims

    if self.test_mode:
      LOGGER.setLevel(DEBUG)
    else:
      LOGGER.setLevel(INFO)

  def make_sample_sheet(self, libs, fname):

    '''Given a dict of adapters (library code keys, adapter sequence
    values) and the source fastq filename, write out a sample sheet
    suitable for use with demuxIllumina.'''

    (_limscodes, flowcell, flowlane, flowpair) = \
        parse_incoming_fastq_name(fname)

    sample_fn = "sampleSheet_%d.txt" % (os.getpid(),)
    fdesc = open(sample_fn, 'w')

    for code in libs:
      try:
        lib = Library.objects.get(code=code)
      except Library.DoesNotExist, _err:
        raise StandardError("Unable to find library code %s in the database."
                            % code)
      if not lib.adapter:
        raise StandardError("Library has no adapter associated: %s." % code)
      fname = build_incoming_fastq_name(code.lower(), flowcell,
                                        flowlane, flowpair)
      fdesc.write("%s %s\n" % (lib.adapter.sequence, fname))
      if code not in self._demux_files:
        self._demux_files[code] = set()
      self._demux_files[code].add(fname)

    fdesc.close()
    return sample_fn

  def demultiplex(self, codes, fname):
    '''Actually run the demultiplexing, using demuxIllumina.'''

    # look up adapters from database,
    # write sampleSheet file
    LOGGER.debug("Making sample sheet.")
    sheet = self.make_sample_sheet(codes, fname)
    LOGGER.info("Sample sheet created.")
    # invoke demultiplexer
    cmd = [self.demux_prog, '-d', sheet, fname] # demuxIllumina v2.0 and above
    LOGGER.debug("Command for demultiplexing: %s", " ".join(cmd))
    pout  = call_subprocess(cmd, path=self.conf.hostpath)
    fnpat = re.compile(r"tag\s+(\w+):\s+([^\s]+)\s*$")
    fnset = set()
    lostpat = re.compile(r"lost\s+([\/\d]+)\s+reads")
    for line in pout:
      matchobj = fnpat.match(line)
      if matchobj:
        fnset.add(matchobj.group(2))
      else:
        matchobj2 = lostpat.match(line)
        if matchobj2:
          LOGGER.info("lost %s reads", matchobj2.group(1))

    # Calculate MD5 checksums. I'm honestly not clear what we actually
    # use these files for? FIXME
    for fname in fnset:
      md5   = checksum_file(fname)
      md5fn = fname + ".md5"
      with open(md5fn, 'wb') as md5out:
        md5out.write(md5 + "\n")
      set_file_permissions(self.conf.group, fname)
      set_file_permissions(self.conf.group, md5fn)

    # Delete the sample sheet.
    os.unlink(sheet)
    # report results

  def run(self, flowcell, flowlane=None, fcq=None):
    '''The main entry point for the class.'''
    multiplexed = {}

    # get list of new lanes from flow cell
    if fcq is None:
      fcq = FlowCellQuery(flowcell, flowlane, lims=self.lims)

    flowlanes = set()
    if fcq.lims_fc.analysis_status not in self.ready:
      LOGGER.info("flow cell status '%s'", fcq.lims_fc.analysis_status)
      sys.exit("Flow cell analysis status not yet completed.")

    for (lanenum, libset) in fcq.lane_library.items():
      if lanenum not in multiplexed:
        multiplexed[lanenum] = set()
      for lib in libset:
        if fcq.lib_status[lib] in ('new') or not self.db_library_check:

          # Only register lane for demultiplexing if this if lib not
          # in lane.lims_samples()
          if not fcq.lane_demuxed[lanenum]:
            multiplexed[lanenum].add(lib)

          flowlanes.add((fcq.lims_fc.fcid, lanenum))

    if len(flowlanes) == 0:
      LOGGER.info("No ready lanes for flowcell '%s'", flowcell)
      sys.exit("No lanes to process.")

    # We need to set our working directory to something suitable
    # before we start; otherwise we end up demuxing into a home
    # directory or similar.
    pwd = os.getcwd()
    os.chdir(self.conf.incoming)

    downloading = Status.objects.get(code='downloading data')

    # for each lane...
    path = self.conf.incoming
    for (flowcell, flowlane) in flowlanes:

      # Mark our lane(s) as active (note that each library has its own
      # version of this lane).
      for lane in Lane.objects.filter(flowcell=flowcell, flowlane=flowlane):
        lane.status = downloading
        lane.save()

      # retrieve file
      fetcher = FQFileFetcher(destination=path, lims=self.lims,
                              test_mode=self.test_mode)
      fetcher.fetch(flowcell, flowlane)

      if self.test_mode:
        print ("Test Mode: skipping download of %s lane %s to %s"
               % (flowcell, flowlane, path))
        continue

      for fname in fetcher.targets:
        if len(fname) > 0:

          # Check file was retrieved.
          if not os.path.exists(fname):
            LOGGER.error("Can't seem to find expected file '%s'", fname)
          else:
            muxed_libs = multiplexed[flowlane]
            if len(muxed_libs) > 1:

              # Demultiplex file if required.
              LOGGER.info("Demultiplexing file %s for libraries: %s",
                          fname, ", ".join(muxed_libs))
              self.demultiplex(muxed_libs, fname)
              for lib in muxed_libs:
                self.output_files += self._demux_files[lib]
            else:
              LOGGER.info("File does not require demultiplexing: %s", fname)
              self.output_files.append(fname)

    os.chdir(pwd)
    LOGGER.info("Initial FlowCell processing complete.")
    return self.output_files

###############################################################################

class FlowCellQuery(object):

  '''Main query class.'''

  __slots__ = ('verbose', 'lims_fc', 'lib_status', 'lane_library',
               'lane_demuxed', 'lims')

  def __init__(self, flowcell_id, flow_lane=None, verbose=False, lims=None):
    LOGGER.setLevel(INFO)
    self.verbose = verbose

    # Map libcode to status
    self.lib_status   = {}

    # Map lanenum to libraries
    self.lane_library = {}

    # Map lanenum to demux flag (True = demuxed files available)
    self.lane_demuxed = {}

    self.lims_fc    = None
    if lims is None:
      lims = Lims()
    if not lims.running():
      LOGGER.error("Remote LIMS access broken... cannot continue.")
      sys.exit("LIMS not running.")
    self.lims = lims
    self.lims_fc    = self._run_query(flowcell_id, flow_lane)

  def check_lane(self, lane, lims_fc):
    '''Retrieve lane, flowcell and library information from the
    repository database.'''
    if lane.lane not in self.lane_library:
      self.lane_library[lane.lane] = set()
    #  user_code = lane.user_sample_id.lower().replace("-", "").replace(" ", "_")
    codes = demux_code(lane.user_sample_id)
    for code in codes:
      try:
        user_code = LibraryNameMap.objects.get(limsname=code).libname
      except LibraryNameMap.DoesNotExist, _err:
        user_code = code
      try:
        lib = Library.objects.search_by_name(user_code)
      except Library.DoesNotExist, _err:
        print "%s %s %s %s" % (lims_fc.fcid, lane.lane, user_code, 'noLibrary')
        continue

      lib_status = None
      db_lanes = Lane.objects.filter(flowcell=lims_fc.fcid,
                                    flowlane=lane.lane,
                                    library=lib)
      if db_lanes.count() == 1:
        # if db_lane != None:
        db_lane = db_lanes[0]
        db_lib  = db_lane.library
        if lib != db_lib:
          LOGGER.warn("%s %s: LIMS code='%s' but LIB code='%s'",
                      lims_fc.fcid, lane.lane, lib.code, db_lib.code)
        lib_status = 'processed'
        print "%s %s %s %s %s" % (lims_fc.fcid, lane.lane,
                                  db_lib.code, lib_status, db_lane.name)
      else:
        lib_status = 'new'
        print "%s %s %s %s" % (lims_fc.fcid, lane.lane, lib.code, lib_status)

      self.lib_status[lib.code] = lib_status
      self.lane_library[lane.lane].add(lib.code)

      if lib.code.lower() in [ x.lower() for x in lane.lims_demuxed_samples() ]:
        self.lane_demuxed[lane.lane] = True
      else:
        self.lane_demuxed[lane.lane] = False

  def _run_query(self, flowcell_id, flowlane_num):
    '''Main entry point for the class (called via __init__).'''

    # Only run the query once.
    if self.lims_fc:
      return self.lims_fc

    # Only look for active users (otherwise we get lanes from old lab
    # members still working within the building).
    people = User.objects.filter(is_active=True)
    emails = munge_cruk_emails([x.email.lower() for x in people])
    lims_fc = self.lims.load_flowcell(flowcell_id)
    if lims_fc == None:
      sys.exit("Problem retrieving flow cell with ID %s from LIMS."
               % (flowcell_id,))
    if self.verbose:
      lims_fc.dump()

    if lims_fc.fcid != flowcell_id:
      LOGGER.error(
        "Flow cell id in LIMS result %s != %s from command line.  Quitting.",
        lims_fc.fcid, flowcell_id)
      sys.exit("LIMS query retrieved unexpected results.")

    try:
      print ("%s %s %s" % (flowcell_id,
                           lims_fc.analysis_status,
                           lims_fc.finish_date))
    except AttributeError, _err:
      print "%s UNK %s" % (flowcell_id, lims_fc.finish_date)
    for lane in lims_fc.iter_lanes():
      if flowlane_num == None or lane.lane == int(flowlane_num):
        if any([x.lower() in emails for x in lane.user_emails]):
          self.check_lane(lane, lims_fc)

    return lims_fc

