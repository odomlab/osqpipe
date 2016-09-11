#!/bin/env python
#
# $Id$

'''Script to add a lane-associated file to the repository.'''

import sys

import os.path
import datetime

from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG

# New in Django 1.7 and above.
import django

from osqutil.utilities import parse_repository_filename, checksum_file, run_in_communication_host
from osqpipe.models import Filetype, Lane, Lanefile, Facility, ArchiveLocation
from osqpipe.pipeline.laneqc import LaneFastQCReport

###########################################################

def checksum_from_file(fname):
  md5 = None
    
  fnmd5 = fname + '.md5'
  if os.path.exists(fnmd5):
    with open(fnmd5, 'rb') as fh:
      line = fh.readline().rstrip('\n')  
      (md5, tmp) = line.split()    
    
  return md5

def get_lane_for_file(fname):

  lane = None
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

  if lane is None:
    LOGGER.error("Lane not determined! Exiting.")
    sys.exit(1)

  return lane

def date_to_sqldate(self, datestr):
        sqldate = ''
        if len(datestr) == 6:
            sqldate = '20%s-%s-%s' % (datestr[:2], datestr[2:4], datestr[4:6])
        return sqldate

class RepoFileHandler(object):

  '''Class which is almost certainly overkill given the limited
  functionality left in this script, post-refactor.'''

  @staticmethod
  def run(fns, md5files=False, archive=None):
    '''Main entry point for the class.'''

    arc = None
    arc_date = None

    if archive is not None:
      try:
        arc = ArchiveLocation.objects.get(name=archive)
        arc_date = datetime.date.today()
      except ArchiveLocation.DoesNotExist, _err:
        raise SystemExit("No ArchiveLocation with name '%s'" % archive)
      
    for fname in fns:
      # We assume that the file names in the list may correspond to different lanes.
      # Hence, we search lane for each file again.
      lane = get_lane_for_file(fname)
      # Even though fname was already parsed in get_lane_for_file, parse it again as we need the pipeline value
      (code, facility, lanenum, pipeline) = parse_repository_filename(fname)
      
      if md5files:
        chksum = checksum_from_file(fname)
        if chksum is None:        
          chksum = checksum_file(fname)          
      else:
        chksum = checksum_file(fname)
      filetype = Filetype.objects.guess_type(fname)
      basefn = os.path.split(fname)[1]
      fnparts = os.path.splitext(basefn)
      if fnparts[1] == '.gz':
        basefn = fnparts[0]
      LOGGER.debug("basefn: '%s'" % (basefn))    

      lanefile = Lanefile(filename=basefn, checksum=chksum,
                          filetype=filetype, lane=lane,
                          pipeline=pipeline, archive=arc, archive_date=arc_date)
      lanefile.save()
      LOGGER.info("Added %s to repository.", basefn)

  @staticmethod
  def add_qc_files(fnames, program_name):
    # Assume first file is representative of lane for all
    l = get_lane_for_file(fnames[0])

    LOGGER.info("Inserting QC files for lane=%d", l.id)

    with LaneFastQCReport(lane=l, program_name=program_name, workdir='./', move_files=False) as rep:
      rep.output_files = fnames
      rep.insert_into_repository()

  @staticmethod          
  def add_lane_summary(fname):

    if not os.path.exists(fname):
      LOGGER.error("%s is missing or unaccessible", fname)
      sys.exit(1)
    
    lane = get_lane_for_file(fname)

    goodreads = []
    badreads = []
    readlen = 0
    with open(fname, 'rb') as fh:
      for line in fh:
        line = line.rstrip('\n')
        flds = line.split()
        if flds[0] == "empty":
          LOGGER.warning("%s: no data", fname)
          # lane.runnumber = 'unknown'
        keyword = flds[0]
        if len(flds) > 2:
          data = [ float(x) for x in flds[1:] ]
        else:
          data = flds[1]

        if keyword == 'good':
          goodreads.append(data)
          if len(data) > readlen:
            readlen = len(data)
        elif keyword == 'bad':
          badreads.append(data)
        elif keyword in ('runnumber','reads','passedpf',
                           'qualmean','qualstdev','qualmeanpf','qualstdevpf'):

          # Note that we omit flowlane deliberately; it is no longer
          # parsed correctly by summarizeFile (and is not all that
          # desirable to change at this point in the code!). Also,
          # machine is now better identified via the upstream LIMS.
          # This is a bit ugly. Is there a better way using Django?
          vars(lane)[keyword] = data
          LOGGER.debug("LaneInfo: '%s' => '%s'", keyword, flds[1])
          # Extact read length from the legth of qualmean. qualstdev, qualpeanpf, qualstdevpf would do as well.
        else:
          if keyword in ('machine','flowlane'):
            continue
          else:
            LOGGER.error("%s content does not look like from summarizeFile. Keyword %s not recognized!", fname, keyword)
            sys.exit(1)
    if lane.runnumber is None:
      LOGGER.error("No runnumber information parsed from file header.")
      raise(Exception("Problem collecting information from file."))

    lane.seqsamplepf = "\n".join(goodreads[1:100])
    lane.seqsamplebad = "\n".join(badreads[1:100])
    lane.readlength = readlen
    
    lane.save()
    LOGGER.info("Added file summary from %s", fname)

###########################################################

if __name__ == '__main__':

  # Make sure following is run in communication host.
  run_in_communication_host(sys.argv)

  
  import argparse

  PARSER = argparse.ArgumentParser(
    description='Add a list of files associated with'
    + ' a single lane to the repository. If not provided, computes md5sums on run. By default, it is assumed that files are of type lanefile.')

  PARSER.add_argument('files', metavar='<files>', type=str, nargs='*',
                      help='The list of files.')
  PARSER.add_argument('-m', '--md5sum', dest='md5files', action='store_true', help='Assume pre-computed md5sums available in filename.md5.', default=False)
  PARSER.add_argument('-s', '--file_summary', dest='summary_file', type=str, help='File created by summaryFile from a fastq file.', default=None)
  PARSER.add_argument('--archive', dest='archive', type=str, help='Archive (e.g. bamark, ebiark, ark) where the file has been saved. Deafult=none.', default=None)
  PARSER.add_argument('--qcfile', dest='qcfile', action='store_true', help='Files should be inserted as QC files. Files are then moved to repository!', default=False)
  PARSER.add_argument('--program_name', dest='program_name', type=str, help='Program name used for generating qc files. Default=\'fastqc\'.', default='fastqc')

  ARGS = PARSER.parse_args()

  if len(ARGS.files) == 0 and ARGS.summary_file is None:
    PARSER.print_help()
    sys.exit(1)

  LOGGER = configure_logging(level=DEBUG)
  django.setup()

  HND = RepoFileHandler()

  if ARGS.summary_file is not None:
    HND.add_lane_summary(ARGS.summary_file)
  if len(ARGS.files):
    if ARGS.qcfile:      
      HND.add_qc_files(ARGS.files, ARGS.program_name)
    else:
      HND.run(ARGS.files, md5files=ARGS.md5files, archive=ARGS.archive)
