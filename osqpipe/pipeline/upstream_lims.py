#!/usr/bin/env python
#
# $Id$

'''Code to interact with the main Genologics LIMS installation.'''

import sys
import re
import os.path
import logging
import weakref
from datetime import date, timedelta
import xml.etree.ElementTree as ET

import requests

from config import Config
from utilities import munge_cruk_emails
from ..models import LibraryNameMap, User

from setup_logs import configure_logging
LOGGER = configure_logging('lims')

CONFIG = Config()

###############################################################################
def http_download_file(url, local_filename, params=None):
  '''
  Download a remote file URL to a local filename. Returns the local
  filename on success.
  '''
  # See http://stackoverflow.com/a/16696317 for the original source.

  if params is None:
    params = {}

  # Read large files in chunks using stream = True parameter
  req = requests.get(url, stream = True, params=params)
  if req.status_code is not 200:
    LOGGER.error("Failed to download file: %s", str(params))
    raise StandardError("Unable download LIMS file: %s" % str(params))
  with open(local_filename, 'wb') as outfh:
    for chunk in req.iter_content(chunk_size=1024):
      if chunk: # filter out keep-alive new chunks
        outfh.write(chunk)
        outfh.flush()
  return local_filename

def get_lims_run_history(url, since):
  '''
  Query the LIMS for runsWithFilesAttachedSince; return the result of
  running the query response through xml.etree.ElementTree.fromstring.
  '''
  LOGGER.debug("Querying LIMS for runs since %s", since)
  history_url = "%s/runsWithFilesAttachedSince?start=%s" % (url, since)
  req = requests.get(history_url)
  if req.status_code is not 200:
    LOGGER.error("Failed to retrieve runs since date %s.", since)
    raise StandardError("Unable to retrieve LIMS run history: %s" % history_url)
  try:
    root = ET.fromstring(req.text)
  except ET.ParseError, err:
    LOGGER.error("LIMS query returned bad XML.")
    raise StandardError("Bad XML in response from LIMS query: %s" % err)
  return root

def get_lims_run_details(url, run_id):
  '''
  Query the LIMS for fullDetailsOfRun; return the result of
  running the query response through xml.etree.ElementTree.fromstring.
  '''
  LOGGER.info("Querying LIMS for run_id: %s", run_id)
  req = requests.get("%s/fullDetailsOfRun?runId=%s" % (url, run_id))
  if req.status_code is not 200:
    LOGGER.error("Failed to retrieve detail for %s.", run_id)
    raise StandardError("Unable to retrieve LIMS run detail.")
  try:
    root = ET.fromstring(req.text)
  except ET.ParseError, err:
    LOGGER.error("LIMS query returned bad XML.")
    raise StandardError("Bad XML in response from LIMS query: %s" % err)
  return root

def runs_containing_samples(url, libcode):
  '''
  Query the LIMS for runsContainingSamples; return the result of
  running the query response through xml.etree.ElementTree.fromstring.
  '''
  LOGGER.info("Querying LIMS for library: %s", libcode)
  req = requests.get("%s/runsContainingSamples?sampleName=%s" % (url, libcode))
  if req.status_code is not 200:
    LOGGER.error("Failed to retrieve runs for library %s.", libcode)
    raise StandardError("Unable to retrieve run listing from LIMS.")
  try:
    root = ET.fromstring(req.text)
  except ET.ParseError, err:
    LOGGER.error("LIMS query returned bad XML.")
    raise StandardError("Bad XML in response from LIMS query: %s" % err)
  return root

###############################################################################

class LimsLaneFile(object):
  '''
  Thin class used to hold file information from the LIMS.
  '''
  def __init__(self, uri, lims_id, filetype):
    self.uri     = uri
    self.lims_id = lims_id
    self.filetype = filetype

  def __getattribute__(self, key):
    '''Wrap attribute retrieval to seamlessly deref any weakref.ref links.'''
    value = super(LimsLaneFile, self).__getattribute__(key)
    if key == 'lane':
      return value()
    else:
      return value

  def __setattr__(self, key, value):
    '''Wrap attribute setting to seamlessly weakref.ref lane links.'''
    if key == 'lane':
      value = weakref.ref(value)
    super(LimsLaneFile, self).__setattr__(key, value)

class LimsLane(object):
  '''
  Class representing a given flowcell lane as stored within the
  LIMS.
  '''
  lane = None
  sample_process_id = None
  files = None
  samples = None

  def __init__(self, fields, parent):

    # Define our usage of the LIMS API, as for LimsFlowCell
    # itself. Note that the iter_lanes mechanism from LimsFlowCell
    # effectively makes these objects only safe to use as read-only.
    for key in ('lane', 'user_sample_id',
                'user_email', 'genomics_sample_id'):
      self.__dict__[key] = fields.get(key)

    # Dict of file location lists keyed by filetype string
    # (e.g. 'FASTQ').
    self.__dict__['files'] = fields.get('files', {})
    for ftypelist in self.files.values():
      for lfile in ftypelist:
        lfile.__dict__['lane'] = weakref.ref(self)

    # Dict of file location lists keyed by sample name (which should,
    # in theory, be library code).
    self.__dict__['samples'] = fields.get('samples', {})
    for liblist in self.samples.values():
      for lfile in liblist:
        lfile.__dict__['lane'] = weakref.ref(self)

    self.__dict__['flowcell'] = weakref.ref(parent)

  def __getattribute__(self, key):
    '''
    Wrap attribute retrieval to seamlessly deref any weakref.ref links.
    '''
    value = super(LimsLane, self).__getattribute__(key)
    if key == 'flowcell':
      return value()
    else:
      return value

  def __setattr__(self, key, value):
    '''
    Wrap attribute setting to seamlessly weakref.ref flowcell links.
    '''
    if key == 'flowcell':
      value = weakref.ref(value)
    super(LimsLane, self).__setattr__(key, value)

  def __cmp__(self, other):
    if other == None:
      return 1
    return cmp(self.__dict__['lane'], other.__dict__['lane'])

  def __str__(self):
    return "%d" % (self.__dict__['lane'],)

  def build_summary_url(self):
    '''
    Returns the predicted summary URL for this lane in the upstream
    LIMS web UI.
    '''
    # Introspect to find cached FASTQ file location URIs.
    files = self.files.get('FASTQ', [])

    # Generate a directory URL from the file url
    if len(files) > 0:
      file_url = files[0].uri
      url      = re.sub(r'primary/[^\/]+', '', file_url)
      LOGGER.debug("Generated Summary URL: %s", url)
    else:
      url = ""

    return url

  def lims_files(self, ftype='FASTQ'):
    '''
    List the files for this lane (non-demultiplexed) for a given
    file type.
    '''
    return self.files.get(ftype, [])

  def lims_samples(self):
    '''
    List the samples which have been multiplexed in this lane.
    '''
    return self.samples.keys()

  def lims_demuxed_files(self, sample):
    '''
    Return the demultiplexed LimsLaneFile objects for a given
    sample name (library code).
    '''
    return self.samples.get(sample, [])

  def lims_demuxed_samples(self):
    '''
    List the samples which have been multiplexed in this lane for
    which demultiplexed fastq files are available.
    '''
    return [ x for x in self.lims_samples()
             if len(self.samples.get(x, [])) > 0 ]

  def file_locations(self, ftype='FASTQ'):
    '''
    Return all the LIMS file URIs for a given file type (default: FASTQ).
    '''
    return [ x.uri for x in self.files.get(ftype, []) ]

  def file_lims_ids(self, ftype='FASTQ'):
    '''
    Return all the LIMS file IDs for a given file type (default: FASTQ).
    '''
    return [ x.lims_id for x in self.files.get(ftype, []) ]

  def dump(self, out=sys.stderr, indent=6):
    '''
    String representation of the object for debugging.
    '''
    out.write("%slane %s:\n" % (" " * indent, self.lane))
    instr = " " * (indent+4)
    for key in self.__dict__.keys():
      if key != 'files':
        out.write("%s%s : %s\n" % (instr, key, self.__dict__[key]))
    out.write("%sFiles:\n" % instr)
    for ftype in self.files:
      out.write("%s  %s : %s\n"
                % (instr, ftype, "; ".join(self.file_locations(ftype))))

class LimsFlowCell(object):
  '''
  Class representing a sequencing flowcell from the LIMS.
  '''
  instrument      = None
  run_number      = None
  analysis_status = None

  def __init__(self, run):
    self.fcid        = run['flowcell_id']
    self.finish_date = run['finish_date']
    self.lanes      = {}
    for lanedata in run['sample_lanes']:
      lane = LimsLane(lanedata, self)
      LOGGER.debug("LIMS: loaded lane %s", lane)
      self.lanes[lane.lane] = lane

    # Rather than dumping everything from the LIMS into our object, we
    # take just what we need; that way we know what we need when the
    # LIMS output changes and can insert appropriate adapter code here.
    for key in ('instrument', 'run_number', 'analysis_status'):
      self.__dict__[key] = run.get(key)

  def __str__(self):
    return self.fcid

  def get_lane(self, lane_id):
    '''
    Retrieve the LimsLane object for a given lane number.
    '''
    found = None
    for lane in self.lanes.values():
      if lane.lane == lane_id:
        found = lane
    if found == None:
      LOGGER.info("FC: lane '%s' not present in '%s'", lane_id, self.fcid)
    return found

  def get_sample_lane(self, lane_id, sample_id):
    '''
    Retrieve the LimsLane object for a given lane number and sample
    ID (library code).
    '''
    sample_id = sample_id.lower()
    found = None
    for lane in self.lanes.values():
      try:
        lanelib = LibraryNameMap.objects.get(limsname=lane.user_sample_id).libname
      except LibraryNameMap.DoesNotExist, _err:
        lanelib = lane.user_sample_id
      lanelib = lanelib.lower()

      if lane.lane == lane_id and lanelib == sample_id:
        found = lane
      else:
        LOGGER.debug("checking for '%s' in '%s'", sample_id, lane.user_sample_id)
        flds = [ x.strip().lower() for x in lane.user_sample_id.split(',') ]
        if sample_id in flds:
          found = lane
    if found == None:
      LOGGER.info("FC: lane '%s/%s' not present in '%s'",
                   lane_id, sample_id, self.fcid)
    return found

  def iter_lanes(self):
    '''
    Iterate over all the lanes for this flowcell object.
    '''
    lane_ids = [ (x.lane, x.user_sample_id, x.lane) for x in self.lanes.values() ]
    lane_ids.sort()
    for lane in lane_ids:
      yield self.lanes[lane[2]]
    return

  def dump(self, out=sys.stderr, indent=4):
    '''
    String representation of the object for debugging.
    '''
    out.write("%sflowcell %s:\n" % (" " * indent, self.fcid))
    indent += 4
    instr = " " * indent
    for key in self.__dict__.keys():
      if key != 'lanes':
        out.write("%s%s : %s\n" % (instr, key, self.__dict__[key]))
    for lane in self.iter_lanes():
      lane.dump(out, indent-2)

class Lims(object):
  '''
  Main class representing the Genologics LIMS itself.
  '''
  def __init__(self, uri=None):

    sys.setrecursionlimit(1000)

    if uri == None:
      uri = CONFIG.lims_rest_uri
    self.uri = uri
    LOGGER.debug("LIMS: Pointing at REST API '%s'", uri)

    # Cache for downloaded data (see load_run).
    self._run_details = {}

    # Retrieve the boilerplate help page as a test.
    req = requests.get(self.uri)
    if req.status_code is not 200:
      LOGGER.error("LIMS: REST API is not responding.")
      raise StandardError("Unable to query LIMS REST API.")

  def __str__(self):
    return "LIMS: %s" % self.uri

  def running(self):
    '''
    Test to make sure the LIMS is running and we can connect to it.
    '''
    # Retrieve the boilerplate help page as a test.
    req = requests.get(self.uri)
    if req.status_code is 200:
      return True
    else:
      LOGGER.warning("LIMS: REST API is not responding.")
      return False

  def load_flowcell(self, flowcell):
    '''
    Retrieve flowcell data from the LIMS. In time it is hoped to
    move away from using this and towards using load_run, below.
    '''
    # This currently needs to find the runID for a given flowcell and
    # then call load_run().

    # FIXME we would ideally deprecate this method at some stage in
    # favour of using load_run directly. Alternatively, the upstream
    # LIMS could implement an API search call by flowcell ID.

    runid = self.run_id_from_flowcell(flowcell)
    return self.load_run(runid)

  def run_id_from_flowcell(self, flowcell):
    '''
    Runs a set of queries, up to a currently hard-coded time limit
    of about one month in the past, to map a flowcell ID onto a more
    canonical run ID.
    '''
    runid = None
    fib   = [1, 1]

    while fib[1] < 35: # Last fibonacci no. will actually be 34.

      # Retrieve the run listing for the last fib[1] days. Note that
      # we may want to consider using completeRunsByFinishDate here
      # instead? Although runsWithFilesAttachedSince will give a more
      # reliable runFolder value.
      since = date.today() - timedelta( fib[1] )
      root  = get_lims_run_history(self.uri, since)

      # XPath query to retrieve the run by flowcell ID.
      run_elem = root.find(".//flowcell[flowcellId='%s']/.." % flowcell)
      if run_elem is not None:
        folder = run_elem.find('./runFolder')

        # This is only a problem with runsWithFilesAttachedSince.
        if folder is None:
          raise StandardError(
            "LIMS output unexpectedly lacks runFolder element.")

        runid = folder.text
        break

      fib = [ fib[1], sum(fib) ]

    if runid is None:
      LOGGER.error("No run ID found corresponding to flowcell: %s", flowcell)
      raise StandardError("Unable to retrieve flowcell Run ID.")

    return runid

  def load_run(self, run_id, requery=False):
    '''
    Retrieve the full details for a given run ID. This method will
    cache details and avoid requerying unless forced to do so, as this
    can be a rather slow and presumably expensive query on the server.
    '''
    # Check there isn't already some cached data for this run_id.
    if run_id in self._run_details and not requery:
      LOGGER.debug("Returning cached object for run_id: %s", run_id)
      return LimsFlowCell(self._run_details[run_id])

    root = get_lims_run_details(self.uri, run_id)

    # Now we pull out various bits of information and cache them for
    # later.
    run_elem = root.find('.//run')
    run = {
      'instrument'  : run_elem.find('./instrument').text,
      'run_number'  : run_elem.find('./runFolder').text,
      'flowcell_id' : run_elem.find('./flowcell/flowcellId').text,
      'finish_date' : run_elem.find('./finishDate').text,
      }
    lanes = []
    for lanelib in run_elem.findall('./flowcell/library'):
      lanenum = lanelib.find('./lane').text
      LOGGER.debug("Processing LIMS data from lane %s", lanenum)

      # FIXME at some point we should make user_sample_id a proper list,
      # rather than fudging with a comma-separated string.
      lanedict = {
        'lane'             : int(lanenum),
        'user_sample_id'   : ",".join([ x.text for x in
                                        lanelib.findall('./sample/name')]),
        'genomics_sample_id' : lanelib.find('./slxId').text,
        }

      # There's some inconsistency in the LIMS output regarding where
      # the project element gets attached.
      try:
        lanedict['user_email'] = \
            lanelib.find('./project/owner/email').text
      except AttributeError, _err:
        try:
          lanedict['user_email'] = \
              lanelib.find('./sample/project/owner/email').text
        except AttributeError, _err:
          lanedict['user_email'] = None # Some lanes just don't have an email.

      # Store Lane-level data file info.
      filedict = {}
      ftype_map = {
        re.compile(r'FASTQC Lane Report$')  : 'LANE_FASTQC',
        re.compile(r'FASTQ MD5 Checksums$') : 'FASTQ_MD5',
        re.compile(r'MGA Lane Report$')     : 'LANE_MGA',
        re.compile(r'Read \d+ FASTQ$')      : 'FASTQ',
      }
      for file_elem in lanelib.findall('./file'):
        fdesig = file_elem.find('./artifactName').text
        for ( elem_re, ftype ) in ftype_map.iteritems():
          if elem_re.search(fdesig) is not None:
            LOGGER.debug('Classifying LIMS file artifact %s as %s',
                         fdesig, ftype)
            ftype_list = filedict.setdefault(ftype, [])
            newfile = LimsLaneFile(uri=file_elem.find('./url').text,
                                   filetype=ftype,
                                   lims_id=str(file_elem.attrib['fileLimsId']))
            ftype_list.append(newfile)
            break
      lanedict['files'] = filedict

      # Store Sample-level (demultiplexed) data file info (FASTQ only
      # at the moment).
      sampledict = {}
      wanted = re.compile(r'Read \d+ FASTQ$')
      for sample_elem in lanelib.findall('./sample'):
        libcode = sample_elem.find('./name').text
        demux_list = sampledict.setdefault(libcode, [])
        for file_elem in sample_elem.findall('./file'):
          fdesig = file_elem.find('./artifactName').text
          if wanted.search(fdesig) is not None:
            LOGGER.debug('Found demultiplexed FASTQ file %s', fdesig)
            newfile = LimsLaneFile(uri=file_elem.find('./url').text,
                                   filetype='FASTQ',
                                   lims_id=str(file_elem.attrib['fileLimsId']))
            demux_list.append(newfile)

      lanedict['samples'] = sampledict

      lanes.append(lanedict)

    # If a lane-level FASTQ is present for every lane, or a
    # demultiplexed FASTQ present for every sample in every lane, set
    # analysis_status to 'COMPLETE'. This needs to only check the lanes
    # we're actually interested in.
    people = User.objects.filter(is_active=True)
    emails = munge_cruk_emails([x.email.lower() for x in people])
    ourlanes = [ ldict for ldict in lanes
                        if ldict['user_email'] is not None
                       and ldict['user_email'] in emails ]
    demux_fastq = []
    lane_fastq  = []
    for ldict in ourlanes:
      samples_demuxed = []
      for sdict in ldict['samples'].values():

        # Record presence of demuxed FASTQ for each sample in the lane.
        samples_demuxed.append(any([ lfile.filetype == 'FASTQ'
                                     for sdict in ldict['samples'].values()
                                     for lfile in sdict ]))

      # Record whether all samples have a demuxed FASTQ in this lane.
      demux_fastq.append(all(samples_demuxed))

      # Record whether this lane has a lane-level FASTQ file.
      lane_fastq.append(any([ lfile.filetype == 'FASTQ'
                              for fdict in ldict['files'].values()
                              for lfile in fdict ]))

    if all(lane_fastq) or all(demux_fastq):
      run['analysis_status'] = 'COMPLETE'
    else:
      run['analysis_status'] = 'INCOMPLETE'

    run['sample_lanes'] = lanes
    self._run_details[run_id] = run

    return LimsFlowCell(run)

  def get_file_by_id(self, file_lims_id, local_filename):
    '''
    Download a LIMS file (using its LIMS ID) to a local
    filename. Returns the local filename on success.
    '''
    download_uri = "%s/downloadFile" % self.uri
    return http_download_file(download_uri,
                            local_filename,
                            {'fileLimsId':file_lims_id})

###############################################################################

if __name__ == "__main__":
  LOGGER.setLevel(logging.DEBUG)
  if len(sys.argv) < 2:
    print "usage: %s <flowcell>" % (os.path.basename(sys.argv[0]),)
  else:
    LIMS = Lims()
    FLOW = LIMS.load_flowcell(sys.argv[1])
    if FLOW != None:
      FLOW.dump(sys.stderr)
