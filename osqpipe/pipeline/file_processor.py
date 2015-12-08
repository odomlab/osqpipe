#!/usr/bin/env python
#
# $Id$

'''Script to align fastq files against their genome as specified in
the repository, using bwa.'''

import sys
import os
import os.path
import re
from datetime import date
from shutil import move
from tempfile import mkstemp

from utilities import parse_incoming_fastq_name, call_subprocess, \
    checksum_file, parse_repository_filename, rezip_file, \
    set_file_permissions, get_filename_libcode
from config import Config
from ..models import Filetype, Library, Lane, Lanefile, Facility, \
    Status, LibraryNameMap, Machine
from fastq_aligner import FastqBwaAligner, FastqTophatAligner
from upstream_lims import Lims
from fetch_mga import fetch_mga
from laneqc import LaneFastQCReport

from django.db import transaction

from setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging('file_processor')

CONFIG = Config()

###############################################################################

def parse_sanger_metadata(fnamemeta):
  '''
  Parse information from our custom metadata file format as
  downloaded from the Sanger sequencing service.
  '''
  code = None
  metadata = {}
  lnumber = -1
  attribute = None
  value = None
  pattern = re.compile(r'(do\d+).*', re.I)
  LOGGER.debug("Parsing metadata file %s ...", fnamemeta)
  for line in open(fnamemeta, "r"):
    line = line.rstrip('\n')
    lnumber = lnumber + 1
    if lnumber == 5:
      lnumber = 1
    if lnumber == 1:
      (_key, val) = line.split(": ", 2)
      attribute = val
    if lnumber == 2:
      (_key, val) = line.split(": ", 2)
      value = val
      if metadata.get(attribute, False):
        LOGGER.warning("Repeated attribute %s in file %s.  Quitting.",
                       attribute, fnamemeta)
        sys.exit("File content error.")
      else:
        metadata[attribute] = value
        code_match = pattern.match(value)
        if code_match:
          tmpcode = code_match.group(1).lower()
          if code is None or code is not None and tmpcode == code:
            code = tmpcode
          else:
            LOGGER.warning(
              "Multiple DO numbers (%s and %s) found in file %s.  Quitting.",
              code, tmpcode, fnamemeta)
            sys.exit("File content error. Multiple do numbers found!")
  return (metadata, code)

def get_fastq_readlength(fobj):
  '''
  Figure out the read length from the passed models.LaneFile object.
  '''
  if fobj.filetype.gzip:
    from gzip import GzipFile
    gen = lambda x: GzipFile(x)
  else:
    gen = lambda x: open(x)

  # Currently just assumes that the second line is the first read, and
  # that it is representative.
  rlen = None
  with gen(fobj.repository_file_path) as reader:
    for _num in range(2):
      line = reader.next()
    rlen = len(line.rstrip('\n'))

  return rlen

###############################################################################

class GenericFileProcessor(object):

  '''Base class for all FileProcessor classes.'''

  def __init__(self, fname, fname2, paired, facility,
               options=None, notes=None, test_mode=False,
               libcode=None, flowcell=None, flowlane=None):
    self.test_mode = test_mode
    self.incoming = fname
    self.paired = paired
    self.fname2 = fname2
    (self.basename, self.extension) = os.path.splitext(fname)
    self.files = [fname]
    if paired:
      self.files.append(fname2)
    self.tempfiles = []
    self.outfiles = []
    self.notes = notes
    self.library = None

    if options is None:
      options = {}
    self.options = options

    try:
      (self.libcode, self.flowcell, self.flowlane, _flowpair)\
        = parse_incoming_fastq_name(self.basename, ext='')
      self.flowlane = int(self.flowlane)
    except StandardError, _err:
      if libcode is None or flowcell is None or flowlane is None or flowlane == 0:
        raise StandardError(
          "Problem identifying libcode/flowcell/flowlane:\n%s; consider supplying lane identifying metadata manually." % _err)
      self.libcode  = libcode
      self.flowcell = flowcell
      self.flowlane = flowlane

    LOGGER.debug("Identified %s: %s (lane %d) from fastq filename.",
                 self.libcode, self.flowcell, self.flowlane)

    try:
      self.library = Library.objects.get(code=self.libcode)
    except Library.DoesNotExist, _err:
      raise StandardError("Failed to find library for code '%s'."
                          % self.libcode)

    self.facility = facility
    self.make_lane(facility)

    if self.library.paired != self.paired:
      LOGGER.warning("Library paired/single-end annotation disagrees"
                     + " with number of files being processed.")

  def lane_is_empty(self):
    '''
    Tests whether the current lane contains an empty file.
    '''
    if self.test_mode:
      return False
    fname = self.files[0]
    fdesc = open(fname)
    line = fdesc.readline()
    fdesc.close()
    return line == '' # standard test for EOF: readline returns empty string

  def strip_bar_code(self, fname):
    '''
    Trim off the first few bases of the reads in the file, based on
    the barcode associated with its library.
    '''
    if self.library.barcode != None:
      hlen = len(self.library.barcode)
      tmpnam = fname+".tmp"
      cmd = ('trimFastq', '-h', str(hlen), fname, tmpnam)
      LOGGER.debug(" ".join(cmd))
      if not self.test_mode:
        call_subprocess(cmd, path=CONFIG.hostpath)
      LOGGER.debug("mv %s %s", tmpnam, fname)
      if not self.test_mode:
        move(tmpnam, fname)

  def raw2fastq(self, flag=None):
    '''
    Convert raw qseq data to fastq format.
    '''
    # convert to fastq
    newoutfiles = []
    for fname in self.files:
      fastq = os.path.splitext(fname)[0] + '.fq'
      cmd = [ 'export2fastq' ]
      if flag is not None:
        cmd = cmd + [flag]
      # decide whether to filter (we tend to prefer to do so in cases
      # where projects disagree).
      if any([ x.filtered for x in self.library.projects.all() ]):
        cmd = cmd + ['-f']
      cmd = cmd + [fname, fastq]
      LOGGER.debug(cmd)
      if not self.test_mode:
        pout = call_subprocess(cmd, path=CONFIG.hostpath)
        line = pout.readline()
        (total, kept) = line.split()
        self.reads = int(total)
        self.passedpf = int(kept)
        pout.close()
        # PAIRED_END: only remembers last file's PE count
      newoutfiles.append(fastq)
      self.tempfiles.append(fname)
    self.outfiles = newoutfiles

  def trim_fastq(self):
    '''
    Trim the fastq file as specified by the trimhead and trimtail
    options.
    '''
    tnames = []
    for fname in self.outfiles:
      (base, ext) = os.path.splitext(fname)
      headtrim = int(self.options.get('trimhead', 0))
      tailtrim = int(self.options.get('trimtail', 0))
      trimname = "%s_h%d_t%d%s" % (base, headtrim, tailtrim, ext)
      cmd = ('trimFastq', '-h', str(headtrim),
                          '-t', str(tailtrim), fname, trimname)
      LOGGER.debug(" ".join(cmd))
      if not self.test_mode:
        call_subprocess(cmd, path=CONFIG.hostpath)
      tnames.append(trimname)
      self.tempfiles.append(trimname)
    return tnames

  def convert_solexa2phred(self):
    '''
    Convert the Solexa QC scoring system to the Sanger phred score.
    '''
    for fname in self.files:
      if os.path.splitext(fname)[1] == ".fq":
        LOGGER.info("Converting quality values to Sanger format ...")
        tmpname = fname + ".orig"
        LOGGER.debug("mv %s %s", fname, tmpname)
        if not self.test_mode:
          move(fname, tmpname)
        cmd = ('solexa2phred', tmpname, fname)
        LOGGER.debug(" ".join(cmd))
        if not self.test_mode:
          call_subprocess(cmd, path=CONFIG.hostpath)
        self.tempfiles.append(tmpname)

  def make_lane(self, facility):
    '''
    Retrieve or create a new repository lane to hold our data.
    '''
    facobj    = Facility.objects.get(code=facility)

    # Status must be anything *except* the CONFIG.core_ready_status
    # value, so this lane being processed is not picked up again by
    # our cron job.
    newstatus = Status.objects.get(code='processing', authority=None)

    # This needs to handle pre-existing lanes without data (and
    # increment lanenum appropriately).
    try:
      self.lane = Lane.objects.get(facility=facobj,
                                   library=self.library,
                                   flowcell=self.flowcell,
                                   flowlane=self.flowlane)

      # FIXME we don't want a long-range transaction to manage it, but
      # it would be nice if subsequent failure in this pipeline reset
      # the status back to whatever we started with.
      self.lane.status = newstatus
      self.lane.save()

    except Lane.DoesNotExist, _err:

      # Create new lane object.
      self.lane = Lane(facility = facobj,
                       library  = self.library,
                       flowcell = self.flowcell,
                       flowlane = int(self.flowlane),
                       lanenum  = Lane.objects.next_lane_number(self.library),
                       status   = newstatus)

    # Update object with other attrs.
    self.lane.paired = self.paired
    self.lane.notes  = self.notes

  def rename_files(self, facility):
    '''Rename LIMS files according to our current naming scheme.'''
    newnames = []
    for fname in self.files:
      try:
        (_libcode, _flowcell, _flowlane, flowpair)\
            = parse_incoming_fastq_name(fname)
        if self.paired:
          newfn = "%s_%s%02dp%s%s" % (self.library.filename_tag,
                                      facility,
                                      self.lane.lanenum,
                                      flowpair,
                                      self.extension)
        else:
          newfn = "%s_%s%02d%s" % (self.library.filename_tag,
                                   facility,
                                   self.lane.lanenum,
                                   self.extension)

        # Just replace spaces for now as e.g. UCSC upload fails in
        # these cases. Also forward slashes. And parentheses/semicolons.
        sanity_re = re.compile(r'([ \/\(\);]+)')
        newfn     = sanity_re.sub('_', newfn)

        LOGGER.debug("mv %s %s", fname, newfn)
        if not self.test_mode:
          move(fname, newfn)
        newnames.append(newfn)

      except StandardError, _err:
        newnames.append(fname)
    self.files = newnames

  def collect_lane_info(self, flag=None):
    '''
    Retrieve some lane metadata directly from the fastq file.
    '''
    cmd = ['summarizeFile']
    if flag is not None:
      cmd = cmd + [ flag ]
    cmd = cmd + [ self.files[0] ]
    LOGGER.debug(" ".join(cmd))
    if self.test_mode:
      return
    pout = call_subprocess(cmd, path=CONFIG.hostpath)
    goodreads = []
    badreads = []
    for line in pout:
      flds = line.split()
      if flds[0] == "empty":
        LOGGER.warning("%s: no data", self.files[0])
        self.lane.runnumber = 'unknown'
      keyword = flds[0]
      if len(flds) > 2:
        data = [ float(x) for x in flds[1:] ]
      else:
        data = flds[1]

      if keyword == 'good':
        goodreads.append(data)
      elif keyword == 'bad':
        badreads.append(data)
      elif keyword in ('runnumber','reads','passedpf',
                       'qualmean','qualstdev','qualmeanpf','qualstdevpf'):

        # Note that we omit flowlane deliberately; it is no longer
        # parsed correctly by summarizeFile (and is not all that
        # desirable to change at this point in the code!). Also,
        # machine is now better identified via the upstream LIMS.

        # This is a bit ugly. Is there a better way using Django?
        vars(self.lane)[keyword] = data
        LOGGER.debug("LaneInfo: '%s' => '%s'", keyword, flds[1])

    pout.close()
    if self.lane.runnumber is None:
      LOGGER.error("No runnumber information parsed from file header.")
      raise(Exception("Problem collecting information from file."))
#    if 'runnumber' in vars(self.lane):
#      chars_in_num = sum([ not x.isdigit() for x in self.lane.runnumber ])
#      if chars_in_num > 0:
#        # must be an old-style lane, with flow cell instead of run number
#        self.lane.flowcell = self.lane.runnumber
#        self.lane.runnumber = None
    self.lane.seqsamplepf = "\n".join(goodreads[1:100])
    self.lane.seqsamplebad = "\n".join(badreads[1:100])

  def collect_lims_info(self):
    '''
    Retrieve lane metadata from our upstream LIMS.
    '''
    lims_fc = None
    if self.facility != 'SAN':
      lims = Lims()
      lims_fc = lims.load_flowcell(self.flowcell)
    if lims_fc == None: # Typically Sanger pipeline
      if self.lane.rundate is None:
        self.lane.rundate = date(2008, 1, 1)
      return
    self.lane.rundate = lims_fc.finish_date

    # This will already raise an error if machine not found, no
    # further try-catch required. Also, coding it like this allowed me
    # to figure out why it was failing (Machine was not being
    # imported). Don't just catch all exceptions, be specific.
    self.lane.machine = Machine.objects.get(code__iexact=str(lims_fc.instrument))
    lims_lane = lims_fc.get_sample_lane(self.flowlane, self.libcode)
    if lims_lane != None:
      self.lane.usersampleid = lims_lane.user_sample_id
      self.lane.genomicssampleid = lims_lane.genomics_sample_id
      self.lane.summaryurl = lims_lane.build_summary_url()

  def collect_info_empty(self):
    '''
    Dummy method to fill in as much metadata we can for lanes without
    a functional fastq fil.
    '''
    lims_fc = None
    if self.facility != 'SAN':
      lims = Lims()
      lims_fc = lims.load_flowcell(self.flowcell)
    if lims_fc == None: # Typically Sanger pipeline
      if self.lane.rundate is None:
        self.lane.rundate = date(2008, 1, 1)
      if self.lane.machine is None:
        self.lane.machine = Machine.objects.get(code__iexact='unknown')
      self.lane.reads    = 0
      self.lane.passedpf = 0
      self.lane.seqsamplepf  = ''
      self.lane.seqsamplebad = ''
      return
    lims_lane = lims_fc.get_sample_lane(self.flowlane, self.libcode)
    self.lane.usersampleid = lims_lane.user_sample_id
    self.lane.genomicssampleid = lims_lane.genomics_sample_id
    self.lane.rundate = lims_fc.finish_date
    self.lane.machine = Machine.objects.get(code__iexact=str(lims_fc.instrument))
    self.lane.runnumber = lims_fc.run_number
    self.lane.flowlane = self.flowlane
    self.lane.seqsamplepf = ''
    self.lane.seqsamplebad = ''
    self.lane.reads = '0'
    self.lane.passedpf = '0'
    self.lane.mapped = '0'
    self.lane.qualmean = [0.0]
    self.lane.qualstdev = [0.0]
    self.lane.qualmeanpf = [0.0]
    self.lane.qualstdevpf = [0.0]
    if lims_lane != None:
      self.lane.usersampleid = lims_lane.user_sample_id
      self.lane.genomicssampleid = lims_lane.genomics_sample_id
      self.lane.summaryurl = lims_lane.build_summary_url()

  def kick_off_alignment(self, tfiles):
    '''
    Start the alignment of our fastq files against the specified
    genome.
    '''
    if 'noalign' in self.options:
      LOGGER.info("Skipping alignment.")
      return

    nocc = None

    # self.library.factor could be 'None' here...
    if str(self.library.factor) in CONFIG.reallocation_factors:
      nocc = CONFIG.nonuniquereads

    # FIXME we ought to remove 'incoming' as hardcoded here.
    repo_incoming = os.path.join(CONFIG.repositorydir, 'incoming')

    # If RNA-Seq, align using tophat2. If not, use our default bwa.
    if self.library.libtype.code == 'rnaseq':
      aligner = FastqTophatAligner(test_mode=self.test_mode,
                                   samplename=self.library.sample.name,
                                   finaldir=repo_incoming)
      if nocc is not None:
        LOGGER.warning("Unsupported attempt to run tophat2 with read reallocation.")
    else:
      aligner = FastqBwaAligner(test_mode=self.test_mode,
                                samplename=self.library.sample.name,
                                finaldir=repo_incoming)
    
    aligner.align_standalone(filepaths=tfiles,
                             genome=self.library.genome,
                             nocc=nocc,
                             destnames=self.outfiles)

  def post_process(self):
    '''
    Stub method standing in for the main processing steps to be
    performed after our initial scan through the files for metadata,
    and any preprocessing steps have completed.
    '''
    pass

  def post_process_empty(self):
    '''
    As for post_process but specifically for empty files.
    '''
    self.collect_info_empty()
    self.outfiles = self.files[:]
    return Status.objects.get(code='complete', authority=None)

  def get_mgareport(self):
    '''
    Retrieve the MGA report from the upstream LIMS, and attach it to
    the list of files for this lane.
    '''
    # assuming the first outfile has ".fq" suffix
    fnsuffix = self.outfiles[0]
    fnsuffix = os.path.splitext(fnsuffix)[0] + ".mga"
    try:
      mgafiles = fetch_mga(self.flowcell,
                           self.flowlane,
                           "", fnsuffix)
    except Exception, _err:
      LOGGER.error("Unable to retrieve MGA files. Skipping!")
      return

    for fname in mgafiles:
      ## remove html file. No use in keeping it.
      if os.path.splitext(fname)[1] == ".html":
        os.unlink(fname)
      else:
        LOGGER.debug("MGA file: %s", fname)
        self.outfiles.append(fname)

  def save_files_to_lane(self):
    '''
    Create lanefile objects in the database and move the files into
    the main repository filesystem.
    '''
    for fname in self.outfiles:
      if not self.test_mode:
        chksum = checksum_file(fname)
      else:
        chksum = "thisisatestchecksum"
      filetype = Filetype.objects.get(suffix=os.path.splitext(fname)[1])
      (_label, _fac, _lane, pipeline) = parse_repository_filename(fname)
      fobj = Lanefile(filename=fname, checksum=chksum,
                      filetype=filetype,
                      pipeline=pipeline,
                      lane=self.lane)
      if filetype.gzip:
        if not self.test_mode:
          fname = rezip_file(fname)
      dest = fobj.repository_file_path
      destdir = os.path.dirname(dest)
      if not os.path.exists(destdir):
        os.makedirs(destdir)
      LOGGER.info("mv -i %s %s", fname, dest)
      if not self.test_mode:

        # It's important to save to the database before moving the
        # file; otherwise we might clobber duplicated file names.
        fobj.save()

        # Note that this remains slightly vulnerable to name
        # collisions between file classes, e.g. Lanefile vs. Alnfile
        # (Update: this is no longer the case now that both are
        # subtypes of Datafile).
        move(fname, dest)
        set_file_permissions(CONFIG.group, dest)

        # Get the read length directly from the fastq file.
        if fobj.filetype.code == 'fq':
          try:
            self.lane.readlength = get_fastq_readlength(fobj)
          except IOError, _err:
            LOGGER.warning("Unable to detect read length from fastq file.")

  def clean_up(self):
    '''
    Delete all temporary files.
    '''
    for fname in self.tempfiles:
      LOGGER.info("Deleting temporary file %s", fname)
      if not self.test_mode:
        os.unlink(fname)

  @transaction.atomic
  def save_all_to_database(self):
    '''
    Save all our files and metadata to the repository database. This
    method is wrapped in a db transaction.
    '''
    if not self.test_mode:

      # FIXME note that the lane information here is not as protected
      # by the transaction mechanism as we would like; we should
      # probably create self.lane as a dict and only create/update
      # Lane in the database within the transaction here. In practice,
      # though, changes to Lane are rare.
      if self.lane.id is None:

        # Save the lane data to the database so we can link to it.
        self.lane.save()

    # Create the lanefile objects.
    self.save_files_to_lane()
    self.clean_up()

    if not self.test_mode:

      # Run FastQC (bam, fastq file formats supported). Note that many
      # of our GenericFileProcessor subclasses might fail here, but we
      # just let that happen for now.
      try:
        with LaneFastQCReport(lane=self.lane, path=CONFIG.hostpath) as qcrep:
          qcrep.insert_into_repository()
      except Exception, err:
        LOGGER.warning("FastQC report generation failed: %s", err)

      # One last save to make sure (this will e.g. update the
      # readlength attribute).
      self.lane.save()

###############################################################################

class ChIPQseqFileProc(GenericFileProcessor):
  '''
  Processor for ChIP-Seq qseq files.
  '''
  def post_process(self):
    '''
    Convert to fastq format, strip barcode, trim fastq, align, collect metadata.
    '''
    self.raw2fastq('-q')
    for fname in self.outfiles:
      self.strip_bar_code(fname)
    if 'trimhead' in self.options or 'trimtail' in self.options:
      tfiles = self.trim_fastq()
    else:
      tfiles = self.outfiles[:]
    self.kick_off_alignment(tfiles)
    self.collect_lims_info()
    self.collect_lane_info("-q")
    return Status.objects.get(code='alignment', authority=None)

###############################################################################

class ChIPExportFileProc(ChIPQseqFileProc):
  '''
  Processor for ChIP-Seq export files.
  '''
  def post_process(self):
    '''
    Convert to fastq format, trim fastq, align, collect metadata.
    '''
    self.raw2fastq()
    if 'trimhead' in self.options or 'trimtail' in self.options:
      tfiles = self.trim_fastq()
    else:
      tfiles = self.outfiles[:]
    self.kick_off_alignment(tfiles)
    self.collect_lims_info()
    self.collect_lane_info("-e")
    return Status.objects.get(code='alignment', authority=None)

###############################################################################

class ChIPFastqFileProc(ChIPQseqFileProc):
  '''
  Processor for ChIP-Seq fastq files.
  '''
  def post_process(self):
    '''
    Convert to phred scoring if required, trim fastq, align, collect metadata.
    '''
    self.outfiles = self.files[:]
    if self.options.get('convert', False):
      self.convert_solexa2phred()
    else:
      LOGGER.info("Assuming quality values in Sanger format.")
    if 'trimhead' in self.options or 'trimtail' in self.options:
      tfiles = self.trim_fastq()
    else:
      tfiles = self.outfiles[:]
    for fname in tfiles:
      set_file_permissions(CONFIG.group, fname)
    self.kick_off_alignment(tfiles)
    self.collect_lims_info()
    self.collect_lane_info("-f")
    return Status.objects.get(code='alignment', authority=None)


###############################################################################

class ChIPMaqFileProc(GenericFileProcessor):
  '''
  Processor for ChIP-Seq map files.
  '''
  pass

###############################################################################

class MNaseFastqFileProc(GenericFileProcessor):
  '''
  Processor for MNase-Seq fastq files.
  '''
  def post_process(self):
    '''
    Convert to phred scoring if required, trim fastq, align, collect metadata.
    '''
    self.outfiles = self.files[:]
    if self.options.get('convert', False):
      self.convert_solexa2phred()
    else:
      LOGGER.info("Assuming quality values in Sanger format.")
    if 'trimhead' in self.options or 'trimtail' in self.options:
      tfiles = self.trim_fastq()
    else:
      tfiles = self.outfiles[:]
    for fname in tfiles:
      set_file_permissions(CONFIG.group, fname)
    self.kick_off_alignment(tfiles)
    self.collect_lims_info()
    self.collect_lane_info("-f")
    return Status.objects.get(code='alignment', authority=None)


###############################################################################

class BisulphiteFileProc(ChIPQseqFileProc):
  '''
  Processor for Bisulphite-Seq qseq files.
  '''
  def kick_off_alignment(self, tfiles):  # Using something like Bismark here?
    raise NotImplementedError()

  def post_process(self):
    '''
    Convert fastq format, trim fastq, collect metadata.
    '''
    self.raw2fastq('-q')
    if 'trimhead' in self.options or 'trimtail' in self.options:
      _tfiles = self.trim_fastq()
    else:
      _tfiles = self.outfiles[:]
    self.collect_lims_info()
    self.collect_lane_info("-q")
    return Status.objects.get(code='complete', authority=None)

###############################################################################

class BisulphiteFastqFileProc(ChIPQseqFileProc):
  '''
  Processor for Bisulphite-Seq fastq files.
  '''
  def kick_off_alignment(self, tfiles):  # Using something like Bismark here?
    raise NotImplementedError()

  def post_process(self):
    '''
    Convert to phred QC scoring if required, trim fastq, collect metadata.
    '''
    self.outfiles = self.files[:]
    if self.options.get('convert', False):
      self.convert_solexa2phred()
    else:
      LOGGER.info("Assuming quality values in Sanger format.")
    if 'trimhead' in self.options or 'trimtail' in self.options:
      _tfiles = self.trim_fastq()
    else:
      _tfiles = self.outfiles[:]
    self.collect_lims_info()
    self.collect_lane_info("-f")
    return Status.objects.get(code='complete', authority=None)

###############################################################################

class MiRFastqFileProc(GenericFileProcessor):
  '''
  Processor for miRNA-Seq fastq files.
  '''
  def _derive_fastq_basename(self, fname):
    '''
    Get the base name of the passed in filename, also accounting for
    .gz extensions.
    '''
    parts   = os.path.splitext(fname)

    # On the rare occasion we pass an already gzipped files to these
    # methods we'd like things to work as we expect.
    if parts[1].lower() == '.gz':
      parts = os.path.splitext(parts[0])
      
    return(os.path.basename(parts[0]))
  
  def trim_linkers(self, fname):
    '''
    Remove linker sequences from the fastq file.
    '''
    # With the TruSeq kit we use, this should never be the case.
    assert self.library.linkerset is not None
    
    base    = self._derive_fastq_basename(fname)
    outfile = base + '_scr.fa'

    # Read lengths of 50bp combined with an expected small RNA length
    # of 20bp and a cut-down 3' adapter sequence of 25bp implies we
    # would need a -3p-prefix offset of at least 5.
    cmd = ('reaper',

           # Input formats and adapter sequences:
           '-i', fname, '-geom', 'no-bc',
           '-3pa',  self.library.linkerset.threep[:25], # reaper bug workaround.
           '-tabu', self.library.linkerset.fivep,
           '-basename', base,

           # Match requirement parameters:
           '-3p-global', '14/2/1/5', # len/edit/gap/offset; default is 14/2/1/0

           # Supposedly the false-positive stripping of real data when
           # using the following option is far less of a problem than
           # biasing all our 3' ends to retain any adapter sequence.
           '-3p-head-to-tail', '1', # Remove perfect matches of at least 1nt.

           # Other, less potentially controversial QC parameters:
           '-nnn-check', '3/5', '-clean-length', '16',
           '--bcq-early',

           # Output format:
           '--nozip', '--noqc',
           '-format-clean', '>i%I_l%L_t%T%n%C%n')
    
    LOGGER.info("Running reaper on %s", fname)
    LOGGER.debug(" ".join(cmd))
    if not self.test_mode:
      call_subprocess(cmd, path=CONFIG.hostpath)
    os.rename("%s.lane.clean" % base, outfile)
    self.tempfiles.append("%s.lint" % base)
    return outfile

  def cluster_exact_matches(self, fname, base=None):
    '''
    Run the tally binary to count exact matches and store reads with
    counts in the output fasta file.
    '''
    if base is None:
      base = self._derive_fastq_basename(fname)
      
    clust_fn = base + ".fa"

    # tally seems to want to compress the outputs even when we don't
    # want it to. FIXME at some point (note that using '-' to indicate
    # STDOUT is also an undocumented tally feature, as far as I can
    # tell). FIXME we may want to remove %L and %T; see whether
    # Nenad's analysis pipeline uses either of them downstream.
    cmd = ('tally -o - -tri 35 -format ">smRNA_%I:count_%C length=%L;trinuc=%T%n%R%n" --fasta-in'
           + ' -i %s | gzip -dc > %s' % (fname, clust_fn))
    LOGGER.info("Running tally on %s", fname)
    LOGGER.debug(" ".join(cmd))
    if not self.test_mode:
      call_subprocess(cmd, path=CONFIG.hostpath, shell=True)
    return clust_fn

  def post_process(self):
    '''
    Convert to phred QC scoring if required, collect metadata, strip
    barcode, trim linkers, convert to fasta, cluster exact matches.
    '''
    if self.options.get('convert', False):
      self.convert_solexa2phred()
    else:
      LOGGER.info("Assuming quality values in Sanger format.")
    self.collect_lane_info()
    self.collect_lims_info()
    self.outfiles = []
    for fname in self.files:
      self.outfiles.append(fname)
      self.strip_bar_code(fname)
      screened_fn = self.trim_linkers(fname)
      cluster_fn = self.cluster_exact_matches(screened_fn,
                                              self._derive_fastq_basename(fname))
      self.kick_off_alignment([cluster_fn])
      self.outfiles.append(cluster_fn)
      self.tempfiles.append(screened_fn)
    return Status.objects.get(code='complete', authority=None)

###############################################################################

class MiRExportFileProc(MiRFastqFileProc):
  '''
  Processor for miRNA-Seq export files.
  '''
  def export2fastq(self, fname, flag=None):
    '''
    Convert input files to fastq format.
    '''
    fq_fn = os.path.splitext(fname)[0] + '.fq'
    cmd = [ 'export2fastq' ]
    if flag is not None:
      cmd = cmd + [flag]
    cmd = cmd + ['-f', fname, fq_fn]
    LOGGER.debug(" ".join(cmd))
    if not self.test_mode:
      call_subprocess(cmd, path=CONFIG.hostpath)
    return fq_fn

  def post_process(self):
    '''
    Collect metadata, convert to fastq format, strip barcode, trim
    linkers, convert to fasta, cluster exact matches.
    '''
    self.collect_lane_info('-e')
    self.collect_lims_info()
    self.outfiles = []
    for fname in self.files:
      fastq_fn = self.export2fastq(fname)
      self.outfiles.append(fastq_fn)
      self.strip_bar_code(fastq_fn)
      screened_fn = self.trim_linkers(fastq_fn)
      cluster_fn = self.cluster_exact_matches(screened_fn,
                                             os.path.splitext(fname)[0])
      self.outfiles.append(cluster_fn)
      self.tempfiles.extend([screened_fn, fname])
    return Status.objects.get(code='complete', authority=None)

###############################################################################

class MiRQseqFileProc(MiRExportFileProc):
  '''
  Processor for miRNA-Seq qseq files.
  '''
  def post_process(self):
    '''
    Collect metadata, convert to fastq format, strip barcode, trim
    linkers, convert to fasta, cluster exact matches.
    '''
    self.collect_lane_info('-q')
    self.collect_lims_info()
    self.outfiles = []
    for fname in self.files:
      fastq_fn = self.export2fastq(fname, '-q')
      self.outfiles.append(fastq_fn)
      self.strip_bar_code(fastq_fn)
      screened_fn = self.trim_linkers(fastq_fn)
      cluster_fn = \
          self.cluster_exact_matches(screened_fn, os.path.splitext(fname)[0])
      self.outfiles.append(cluster_fn)
      self.tempfiles.extend([screened_fn, fname])
    return Status.objects.get(code='complete', authority=None)

###############################################################################

class FileProcessingManager(object):
  '''
  Class of objects used to manage the creation and execution of
  GenericFileProcessor subclasses.
  '''
  __slots__ = ('options', 'facility', 'force_paired_end',
               'libtype2class', 'test_mode')

  def __init__(self, options=None, facility='CRI', force_paired_end=None,
               test_mode=False):
    self.test_mode = test_mode
    if not options:
      options = {}
    self.options  = options
    self.facility = facility
    self.force_paired_end = force_paired_end
    if self.test_mode:
      LOGGER.setLevel(DEBUG)
    else:
      LOGGER.setLevel(INFO)
    self.libtype2class = {
      'chipseq': {'.fq': ChIPFastqFileProc,
                  '.export': ChIPExportFileProc,
                  '.qseq': ChIPQseqFileProc,
                  '.map': ChIPMaqFileProc},
      'ripseq': {'.fq': ChIPFastqFileProc,
                 '.export': ChIPExportFileProc,
                 '.qseq': ChIPQseqFileProc,
                 '.map': ChIPMaqFileProc},
      'chipexo': {'.fq': ChIPFastqFileProc,
                  '.export': ChIPExportFileProc,
                  '.qseq': ChIPQseqFileProc,
                  '.map': ChIPMaqFileProc},
      'chipseq_ssdna_ke': {'.fq': ChIPFastqFileProc,
                           '.export': ChIPExportFileProc,
                           '.qseq': ChIPQseqFileProc,
                           '.map': ChIPMaqFileProc},
      'rnaseq': {'.fq': ChIPFastqFileProc,
                 '.export': ChIPExportFileProc,
                 '.qseq': ChIPQseqFileProc,
                 '.map': ChIPMaqFileProc},
      'exome': {'.fq': ChIPFastqFileProc,
                '.export': ChIPExportFileProc,
                '.qseq': ChIPQseqFileProc,
                '.map': ChIPMaqFileProc},
      'genome': {'.fq': ChIPFastqFileProc,
                 '.export': ChIPExportFileProc,
                 '.qseq': ChIPQseqFileProc,
                 '.map': ChIPMaqFileProc},
      'hic': {'.fq': ChIPFastqFileProc,
              '.export': ChIPExportFileProc,
              '.qseq': ChIPQseqFileProc,
              '.map': ChIPMaqFileProc},
      'mirna': {'.fq': MiRFastqFileProc,
                '.export': MiRExportFileProc,
                '.qseq': MiRQseqFileProc},
      'smrnaseq': {'.fq': MiRFastqFileProc,
                   '.export': MiRExportFileProc,
                   '.qseq': MiRQseqFileProc},
      'ripsmrnaseq': {'.fq': MiRFastqFileProc,
                   '.export': MiRExportFileProc,
                   '.qseq': MiRQseqFileProc},
      'mnaseseq': {'.fq': MNaseFastqFileProc},
      'bisulphite': {'.qseq': BisulphiteFileProc,
                     '.fq': BisulphiteFastqFileProc},
      'bisulph-smrna': {'.fq': MiRFastqFileProc, # Frye lab. Obsolete?
                        '.export': MiRExportFileProc,
                        '.qseq': MiRQseqFileProc}
      }

  def process_sanger_bam(self, fname, metadata, code, library):
    '''
    Convert an incoming bam file from the Sanger sequencing pipeline
    to fastq format for re-alignment using our local reference
    genomes.
    '''
    # Before continuing check that md5 for bam is the same as in
    # metadata file.
    chksum = checksum_file(fname)
    if chksum != metadata['md5']:
      ## Note we cannot assume this will be correct now that we are
      ## postprocessing bam files to remove QC-fail reads.
      LOGGER.warning("Md5sum not same as in metadata file.")

    # As the flowcell id is not known (i.e. is not available easily,
    # improvise flowcell id from run number
    runnumber = fname.split('_')[0]
    lanenumber = fname.split('_')[1][:1]
    # bamname = "%s_SANrun%s_%s.bam" % (code, runnumber, lanenumber)
    bamname = "%s.SANrun%s.s_%s.r.bam" % (code, runnumber, lanenumber)
    if not self.test_mode:
      os.rename(fname, bamname)

    # Extract fastq files
    LOGGER.info("Extract fastq files from bam ...")
    paired = False
    cmd = ('bam2fq', bamname, 'SE')
    if library.paired or self.force_paired_end:
      cmd = ('bam2fq', bamname, 'PE')
      paired = True
    LOGGER.debug(" ".join(cmd))
    if not self.test_mode:
      call_subprocess(cmd, path=CONFIG.hostpath)
    fnametmp1 = "%s_1.fq" % os.path.splitext(bamname)[0]
    fnametmp2 = "%s_2.fq" % os.path.splitext(bamname)[0]
    if not self.test_mode:
      if not os.path.isfile(fnametmp1) and not os.path.isfile(fnametmp2):
        LOGGER.warning(
          "Extraction of fastq files %s %s from %s failed! Quitting.",
          fnametmp1, fnametmp2, bamname)
        sys.exit("Failed to extract fastq files from bam.")

    # Prepare variables to kick off pipeline as usual
    fname = fnametmp1
    fname2 = fnametmp2

    # Following line could in principle be removed as default for
    # convert is already False.
    self.options['convert'] = False
    metaarray = []
    for key in metadata:
      metaarray.append("%s=%s" % (key, metadata[key]))
    notes = ";".join(metaarray)

    return (fname, fname2, notes, bamname, paired)

  @staticmethod
  def retrieve_library(code):
    '''
    Simple wrapper around our database library retrieval code.
    '''
    try:
      library = Library.objects.get(code=code)
    except Library.DoesNotExist, _err:
      LOGGER.warning("No library %s", code)
      sys.exit("No library %s." % (code,))
    return library

  def _run_file_processor(self, fname, fname2, libtype, paired, notes):
    '''
    Create an object of a subclass of GenericFileProcessor and run it
    on our file(s).
    '''
    ext = os.path.splitext(fname)[1]
    libtype2class = self.libtype2class
    if (libtype.code in libtype2class
        and ext in libtype2class[libtype.code]):
      fproc = libtype2class[libtype.code][ext](options=self.options,
                                               fname=fname, fname2=fname2,
                                               paired=paired,
                                               facility=self.facility,
                                               notes=notes,
                                               test_mode=self.test_mode)
    else:
      LOGGER.warning(
        "No class for libtype '%s' and extension '%s'.  Skipping.",
        libtype.code, ext)
      sys.exit("Script does not yet support this library type"
               + " with this filename extension.")

    if fproc.library == None:
      sys.exit("No library code found.")
    fproc.rename_files(self.facility)
    if fproc.lane_is_empty():
      LOGGER.warning("Processing empty file '%s'", fproc.files[0])
      status = fproc.post_process_empty()
    else:
      status = fproc.post_process()
    if self.facility == 'CRI':
      fproc.get_mgareport()

    if status is not None:
      fproc.lane.status = status

    # This is the transaction-wrapped final write to the repository
    # filesystem and database.
    fproc.save_all_to_database()

  def run(self, fns):
    '''The main entry point for this class.'''

    fname  = None
    fname2 = None
    paired = False
    notes  = None
    library = None
    tempfiles = []

    ## check if files are fastq or bam
    ftype = os.path.splitext(fns[0])[1]

    if ftype == '.bam':

      LOGGER.debug("Found bam file. %s", fns[0])
      fname = fns[0]

      if self.facility == 'SAN':

        # For Sanger runs, expect file fname+".meta" in same directory
        # as fname.
        metadata = {}
        fnamemeta = fname + ".meta"

        if os.path.isfile(fname):
          LOGGER.debug("Found metadata file %s.", fnamemeta)
          (metadata, code) = parse_sanger_metadata(fnamemeta)

        # Parse Sanger meta-data file and get library code.
        else:
          LOGGER.warning("Missing '%s' (or not a file).  Quitting.", fnamemeta)
          sys.exit("File access error.")

        library = self.retrieve_library(code)

        # In case of Sanger run, convert bam to fastq files and rename the
        # files according to the library.
        (fname, fname2, notes, bamname, paired) = \
            self.process_sanger_bam(fname, metadata, code, library)

        # Fastq files have been generated! Get rid of Sanger bam file.
        tempfiles.extend([bamname, fnamemeta])

      else:
        LOGGER.error("Unsupported file format: BAM")
        sys.exit("BAM files only supported as part of the Sanger workflow;"
                 + " currently running as %s" % self.facility)

    elif ftype == '.fq' or ftype == '.fastq':
      if len(fns) == 2:
        LOGGER.debug("Found paired reads fastq files: %s, %s", fns[0], fns[1])
        paired = True
        fname = fns[0]
        fname2 = fns[1]
      else:
        LOGGER.debug("Found single read fastq: %s", fns[0])
        fname = fns[0]
        paired = False
        fname2 = None

      # In case of non-Sanger run get the code from fastq filename.
      code = get_filename_libcode(fname)
      try:
        code = LibraryNameMap.objects.get(limsname=code).libname
      except LibraryNameMap.DoesNotExist, _err:
        pass
      library = self.retrieve_library(code)

    else:
      LOGGER.warning(
        "Input file %s of unknown type (based on suffix). Quitting.", fname)
      sys.exit("Unknown input file type error.")

    if not os.path.isfile(fname):
      LOGGER.warning("Missing '%s' (or not a file).  Quitting.", fname)
      sys.exit("File access error.")
    if paired and not os.path.isfile(fname2):
      LOGGER.warning("Missing '%s' (or not a file).  Quitting.", fname)
      sys.exit("File access error.")

    libtype = library.libtype

    # The heavy-lifting part of the code; all downstream processing
    # (e.g., bwa alignment) starts here.
    self._run_file_processor(fname, fname2, libtype, paired, notes)

    if not self.test_mode:

      # Delete temporary files, e.g. Sanger bam and metadata
      for tmpnam in tempfiles:
        os.unlink(tmpnam)


