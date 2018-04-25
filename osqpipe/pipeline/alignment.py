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

'''Script to take a list of bed, bam, bgr and wiggle files, record
them in the database and move the files to the correct part of the
repostory filesystem.'''

import os.path
from shutil import move
import re
import gzip

from django.db import transaction

from osqutil.utilities import is_zipped, parse_repository_filename, \
    checksum_file, rezip_file, set_file_permissions, transfer_file
from ..models import Filetype, Lane, Alignment, Alnfile, Facility, \
    Genome, Program, DataProvenance
from osqutil.samtools import BamToBedConverter
from osqutil.config import Config

from osqutil.progsum import ProgramSummary
from osqutil.setup_logs import configure_logging
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
  checktally = True
  for line in fdesc:
    flds = line.split("\t")
    count = 1

    # Account for tally processing of smallRNA-seq data.
    if checktally:
      tally = flds[3].split(':')
      if len(tally) == 2:
        tally = tally[1].split('_')
        if len(tally) == 2 and tally[0] == 'count':
          count = int(tally[1])
        else:
          checktally = False # these are not tally data.
      else:
        checktally = False # these are not tally data.

    mapped += count
    if int(flds[4]) > 0:
      unique += count
  fdesc.close()
  return (mapped, unique)

###########################################################
class AlignmentHandler(object):

  '''Class designed to manage the insertion of alignment-related files
  into the repository.'''

  __slots__ = ('params', 'prog', 'progvers', 'genome', 'headtrim', 'tailtrim', 'conf')

  def __init__(self, genome, prog, params='', progvers=None, headtrim=0, tailtrim=0):

    # Program and parameters can be a list or scalar. Params elements
    # should always be string; program can be either string or
    # osqpipe.models.Program.
    if all([ type(x) is not list for x in (prog, params) ]):

      # Scalar arguments
      self.prog    = [ prog ]
      self.params  = [ params ]

      # FIXME consider throwing an error here if progvers is already a list.
      self.progvers = [ progvers ]

    elif type(prog) is list:

      # List arguments (params may be the default empty string;
      # progvers may simply by a scalar None)
      self.prog    = prog

      if len(prog) == len(params):
        self.params  = params
      else:
        if params == '': # handle the empty default.
          self.params = [ '' for _x in prog ]
        else:
          raise ValueError("Lengths of prog and params list arguments"
                           + " must match.")
  
      if progvers is None: # handle the empty default.
        self.progvers = [ None for _x in prog ]
      else:
        if len(prog) == len(progvers):
          self.progvers = progvers
        else:
          raise ValueError("Lengths of prog and progvers list arguments"
                           + " must match.")

    else:
      raise TypeError("The params argument cannot be a list if prog is a scalar")

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
    lanelist  = Lane.objects.filter(library__code=code, 
                                    facility__code=facility, 
                                    lanenum=lanenum)
    if lanelist.count() == 0:
      raise ValueError("Could not find lane for '%s'" % (bed))
    elif lanelist.count() > 1:
      raise ValueError(("Found multiple lanes for '%s': "
                       % (bed,)) + ", ".join([x.id for x in lanelist]))
    else:
      lane = lanelist[0]

      aln = self._create_alignment(bed, lane)

    return (aln, lane)

  def _create_alignment(self, bed, lane):

    (mapped, unique) = count_reads(bed)
    gen = Genome.objects.get(code=self.genome)

    # We don't save this yet because we're not currently within a
    # transaction.
    aln = Alignment(lane        = lane,
                    genome      = gen,
                    mapped      = mapped,
                    munique     = unique,
                    total_reads = lane.total_passedpf,
                    headtrim    = self.headtrim,
                    tailtrim    = self.tailtrim)
    return (aln)

  def _find_versioned_program(self, subprog, factor, subvers=None):

    subprog = subprog.strip()

    # If the version has been pre-specified, just use that.
    if subvers is not None:
      subvers = subvers.strip()
      prg = self._retrieve_program_object(program=subprog,
                                          version=subvers)
      return prg

    # If no version specified, look for the program on cluster or althost.
    try:
      althost = self.conf.althost
      assert(althost != '')
    except AttributeError, _err:
      althost = None

    # FIXME come up with a better heuristic than this.
    if subprog in ('reallocateReads', 'samtools') \
          and factor \
          and factor.name in self.conf.reallocation_factors:

      # These programs used on the local server.
      alignerinfo = ProgramSummary(subprog, path=self.conf.hostpath)
        
    else:

      # bwa, maq, gsnap et al. as launched on cluster or remote alignment host.
      if althost is not None and subprog == self.conf.aligner:
        # Using the alternative alignment host.
        alignerinfo = ProgramSummary(subprog, ssh_host=althost,
                                     ssh_user=self.conf.althostuser,
                                     ssh_path=self.conf.althostpath,
                                     ssh_port=self.conf.althostport)
      else:
        # Using the compute cluster as standard.
        alignerinfo = ProgramSummary(subprog, ssh_host=self.conf.cluster,
                                     ssh_user=self.conf.clusteruser,
                                     ssh_path=self.conf.clusterpath,
                                     ssh_port=self.conf.clusterport)
    prg = self._retrieve_program_object(program=alignerinfo.program,
                                        version=alignerinfo.version)
    return prg

  def _retrieve_program_object(self, program, version):
    try:
      prg = Program.objects.get(program=program,
                                version=version,
                                current=True)
    except Program.DoesNotExist, _err:
      raise StandardError("Unable to find current program in database: %s %s"
                          % (program, version))
    return prg

  def _add_data_provenance(self, aln):
    '''
    Update the alignment object with data provenance information.
    '''
    # As program may be made of comma separate list of programs,
    # parse each if these individually and join the results

    for num in range(len(self.prog)):
      subprog = self.prog[num]
      subvers = self.progvers[num]

      # Support the use of string program names or Django Program
      # objects.
      if issubclass(str, type(subprog)):
        prg = self._find_versioned_program(subprog, aln.lane.library.factor, subvers)
      elif issubclass(Program, type(subprog)):
        prg = subprog

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

  @transaction.atomic
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

      # Allow a little latitude when retrieving the file checksum (it
      # may have been stored with or without the trailing .gz).
      try:
        chksum = chksums[fname]
      except KeyError, _err:
        chksum = chksums[basefn]

      # Attempt to save a file record to the database; it's important
      # to do this before moving files, in case there's a naming
      # collision.
      alnfile = Alnfile.objects.create(filename=basefn,
                                       checksum=chksum,
                                       filetype=ftype,
                                       alignment=aln)

      # Move files to permanent locations.
      destname = alnfile.repository_file_path
      # LOGGER.debug("mv %s %s", fname, destname)
      # move(fname, destname)
      # set_file_permissions(self.conf.group, destname)
      LOGGER.debug("mv %s %s@%s:%s", fname, self.conf.user, self.conf.datahost, destname)
      transfer_file(fname, "%s@%s:%s" % (self.conf.user, self.conf.datahost, destname), set_ownership=True)
      os.unlink(fname)

    if final_status is not None:
      aln.lane.status = final_status
      aln.lane.save()

  def add(self, files, final_status=None):

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

    return aln

  def add_bam_to_lane(self, bam, lane, tc1=False, chrom_sizes=None):
    '''
    Generate a bed file from a bam file and add both to the given
    lane. This method is typically used from within an ipython shell to
    handle unusual cases outside the main pipeline. Note that genome
    and data provenance info is passed in via the class attributes
    prog and params.
    '''
    bam_to_bed = BamToBedConverter(tc1=tc1, chrom_sizes=chrom_sizes)
    base       = os.path.splitext(bam)[0]
    
    bedtype = Filetype.objects.get(code='bed')
    bed_fn  = base + bedtype.suffix
    beds    = bam_to_bed.convert(bam, bed_fn)
    chksums = dict( (fname, checksum_file(fname)) for fname in [bam] + beds )

    # First bed file is the main one.
    aln     = self._create_alignment(beds[0], lane)

    if bedtype.gzip:
      bedgz = [ rezip_file(bed) for bed in beds ]

    self._save_to_repository([bam] + bedgz, chksums, aln)
