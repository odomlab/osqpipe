#!/usr/bin/env python
#
# $Id$

'''Script to take a list of bed, bam, bgr and wiggle files, record
them in the database and move the files to the correct part of the
repostory filesystem.'''

import logging
import os.path
from shutil import move
import re
import gzip

from django.db import transaction

from utilities import is_zipped, parse_repository_filename, \
    checksum_file, rezip_file
from ..models import Filetype, Lane, Alignment, Alnfile, Facility, \
    Genome, Program, DataProvenance
from samtools import count_bam_reads
from config import Config

from progsum import ProgramSummary
from setup_logs import configure_logging
LOGGER = configure_logging('alignment')

###########################################################
def count_reads(fname):
  '''
  Count the number of reads in a bed file. This function will handle
  gzipped bed files seamlessly.
  '''
  # FIXME consider using external gzip binary where available.
  if is_zipped(fname):
    fdesc = gzip.open(fname, 'rb')
  else:
    fdesc = open(fname, 'rb')
  mapped = 0
  unique = 0
  for line in fdesc:
    flds = line.split("\t")
    mapped += 1
    if int(flds[4]) > 0:
      unique += 1
  fdesc.close()
  return (mapped, unique)

###########################################################
class AlignmentHandler(object):

  '''Class designed to manage the insertion of alignment-related files
  into the repository.'''

  __slots__ = ('params', 'prog', 'genome', 'headtrim', 'tailtrim', 'conf')

  def __init__(self, genome, prog, params='', headtrim=0, tailtrim=0):
    LOGGER.setLevel(logging.INFO)

    # Program and parameters can be a list, but they must agree on
    # their type. We allow an empty default params string for the sake
    # of a simpler API.
    if (type(prog) != type(params)) and (params != ''):
      raise TypeError("Both prog and params arguments must be of the same"
                      + " type: list, or string.")
    if type(prog) is str:
      self.prog    = [ prog ]
      self.params  = [ params ]
    elif type(prog) is list:
      if len(prog) != len(params):
        if params == '': # handle the empty default.
          self.params = [ '' for _x in prog ]
        else:
          raise ValueError("Lengths of prog and params list arguments"
                           + " must match.")
      self.prog    = prog
      self.params  = params
    else:
      raise TypeError("The prog argument is neither a list nor a string (%s)."
                      % type(prog))

    self.genome  = genome
    self.headtrim = headtrim
    self.tailtrim = tailtrim
    self.conf = Config()

  ## FIXME unused code here?
  @staticmethod
  def calculate_trimming(_fname, name, fq_seq, aln_seq):
    '''Currently unused function?'''
    fq_seq = fq_seq.upper()
    aln_seq = aln_seq.upper()
    if fq_seq == aln_seq:
      return (0, 0)
    pos = fq_seq.find(aln_seq)
    if pos == -1:
      LOGGER.error(
        "%s: aln sequence '%s' does not match fastq sequence '%s'",
        name, aln_seq, fq_seq)
      left = 0
      right = 0
    else:
      left = pos
      right = len(fq_seq) - pos - len(aln_seq)
    return (left, right)

  def aln_from_bedfile(self, bed):

    '''Given a bed file name, retrieve the associated database
    Alignment and Lane objects.'''

    LOGGER.info("Processing bed file: %s", bed)

    # A quick sanity check to try and make sure we don't load a file
    # against the wrong genome. This remains fallible since some
    # genome codes are quite short and might occur in a filename by
    # chance.
    if not re.search(self.genome, bed, re.I): # Note do not merge this re.I change into repackaging branch (it is unnecessary).
      raise ValueError("Genome code not found in bed file name."
                     + " Loading against the wrong genome (%s)?" % self.genome)

    (code, facility, lanenum, _pipeline) = parse_repository_filename(bed)
    lanes  = Lane.objects.filter(library__code=code)
    facobj = Facility.objects.get(code=facility)
    lanelist = [lane for lane in lanes
                if lane.facility == facobj and lane.lanenum == lanenum]
    if len(lanelist) == 0:
      raise ValueError("Could not find lane for '%s'" % (bed))
    elif len(lanelist) > 1:
      raise ValueError(("Found multiple lanes for '%s': "
                       % (bed,)) + ", ".join([x.id for x in lanelist]))
    else:
      lane = lanelist[0]
      (mapped, unique) = count_reads(bed)
      gen = Genome.objects.get(code=self.genome)

      # We don't save this yet because we're not currently within a
      # transaction. Also note that total_reads is not yet set.
      aln = Alignment(lane       = lane,
                      genome     = gen,
                      mapped     = mapped,
                      munique    = unique,
                      headtrim   = self.headtrim,
                      tailtrim   = self.tailtrim)

      # We need to detect the appropriate version. FIXME automatically
      # add new versions to the program table?

    return (aln, lane)

  def _add_data_provenance(self, aln):
    '''
    Update the alignment object with data provenance information.
    '''
    # As program may be made of comma separate list of programs,
    # parse each if these individually and join the results
    try:
      althost = self.conf.althost
      assert(althost != '')
    except AttributeError, _err:
      althost = None
    for num in range(len(self.prog)):
      subprog = self.prog[num]
      subprog.strip()
      if althost is not None and subprog == self.conf.aligner:
        # Using the alternative alignment host.
        alignerinfo = ProgramSummary(subprog, ssh_host=althost,
                                     ssh_user=self.conf.althostuser,
                                     ssh_path=self.conf.althostpath)
      else:
        # Using the compute cluster as standard.
        alignerinfo = ProgramSummary(subprog, ssh_host=self.conf.cluster,
                                     ssh_user=self.conf.clusteruser,
                                     ssh_path=self.conf.clusterpath)
      try:
        prg = Program.objects.get(program=alignerinfo.program,
                                  version=alignerinfo.version,
                                  current=True)
      except Program.DoesNotExist, _err:
        raise StandardError("Unable to find current program in database: %s %s"
                            % (self.prog, alignerinfo.version))

      DataProvenance.objects.create(program      = prg,
                                    parameters   = self.params[num],
                                    rank_index   = num+1,
                                    data_process = aln)

      LOGGER.info("Program=\"%s\"", prg)

  @staticmethod
  def identify_bed_file(files):
    '''From a list of filenames, pick the first bed file we see and
    return it.'''
    bed_re = re.compile('.bed(?:.gz)?$', re.IGNORECASE)
    beds = [ x for x in files if bed_re.search(x) ]
    if len(beds) == 0:
      return None
    else:
      return beds[0]

  @staticmethod
  def identify_bam_file(files):
    '''From a list of filenames, pick the first bam file we see and
    return it.'''
    bam_re = re.compile('.bam$', re.IGNORECASE)
    bams = [ x for x in files if bam_re.search(x) ]
    if len(bams) == 0:
      return None
    else:
      return bams[0]

  @transaction.commit_on_success
  def _save_to_repository(self, files, chksums, aln, final_status=None):
    '''
    Method used to move files into the repository filesystem and
    create Alnfile objects in the database. This is handled within a
    transaction so that problems moving files around don't create
    orphaned database entries.
    '''

    # First, ensure the alignment has been stored in the database.
    aln.save()

    # Propagate the mapped read count up to lane.
    aln.lane.mapped = aln.mapped
    aln.lane.save()

    # Add data provenance from the info stored in self.prog, self.params.
    self._add_data_provenance(aln)

    # Then move the files into position.
    for fname in files:

      ftype = Filetype.objects.guess_type(fname)
      if ftype is None:
        raise ValueError("File type not recognised from database: %s" % fname)

      # Rederive the basefn for recording in the database.
      basefn = os.path.split(fname)[1]
      LOGGER.debug("basefn: '%s'", basefn)
      fnparts = os.path.splitext(basefn)
      if fnparts[1] == self.conf.gzsuffix:
        basefn = fnparts[0]
      LOGGER.debug("basefn: '%s'", basefn)

      # Attempt to save a file record to the database; it's important
      # to do this before moving files, in case there's a naming
      # collision.
      alnfile = Alnfile.objects.create(filename=basefn, checksum=chksums[fname],
                                       filetype=ftype,
                                       alignment=aln)

      # Move files to permanent locations.
      destname = alnfile.repository_file_path
      LOGGER.debug("mv %s %s", fname, destname)
      move(fname, destname)

    if final_status is not None:
      aln.lane.status = final_status
      aln.lane.save()

  def add(self, files, final_status=None, total_reads=None):

    '''
    Process a list of filenames (files must exist on disk). The
    optional final_status argument specifies a
    models.Status object to which the lane should be linked
    upon completion.
    '''

    # We need at least one bed file.
    bed = self.identify_bed_file(files)
    if not bed:
      raise ValueError("Unable to identify any bed files in the input.")

    # Find the appropriate alignment. Note that aln is not yet saved
    # in the database.
    (aln, lane) = self.aln_from_bedfile(bed)

    # Update the alignment with the total number of reads in the
    # supplied bam file. This is mainly to support e.g. smallRNA-seq
    # where the number of reads in the alignment is far smaller than
    # the number in the fastq file.
    if total_reads is None:
      bam = self.identify_bam_file(files)
      if not bam:
        raise ValueError("Unable to identify bam file in input.")
      total_reads = count_bam_reads(bam)
    aln.total_reads = total_reads

    # Do some heavy lifting *outside* of our database transaction, to
    # avoid locking the db for extended periods.
    chksums = dict()
    processed = []
    for fname in files:

      # If the file is uncompressed, don't waste time zipping it prior
      # to checksum.
      chksums[fname] = checksum_file(fname) # also works on zipped file

      ftype = Filetype.objects.guess_type(fname)
      if ftype is None:
        raise ValueError("File type not recognised from database: %s" % fname)
      if ftype.gzip and not is_zipped(fname):
        fname = rezip_file(fname)

      processed.append(fname)

    # All database changes should be handled by the
    # transaction-embedded method below.
    self._save_to_repository(processed, chksums, aln, final_status)

    # In case we want to set lane.status (kept for compatibility;
    # probably defunct though).
    return lane


