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
Script that takes a fastq file and an aligned bam file, creates
a Lane and an Alignment in the repository, filling in the appropriate
values. By default, the fastq file is used only as a source for the
number of reads and is not actually stored in the repository. If the
specified sequencing lane already exists, the script will ignore the
fastq file completely and simply append the bam file to the list of
alignments for that lane.
'''
import os, re
from shutil import move

from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Library, Lane, Alignment, Lanefile,\
    Alnfile, Facility, Genome, Status, Filetype, Program,\
    DataProvenance, Machine
from django.db import transaction
from osqutil.config import Config
from osqutil.utilities import set_file_permissions, is_zipped,\
    rezip_file, unzip_file, checksum_file
from osqpipe.pipeline.laneqc import LaneFastQCReport
from osqpipe.pipeline.alignment import count_reads
from osqutil.samtools import BamToBedConverter
from osqutil.progsum import ProgramSummary

################################################################################

class ExternalDataHandler(object):
  '''
  Class used to add data files to a single lane in the repository. The
  lane in question will be created if it does not already exist.
  '''
  __slots__ = ('config', 'library', 'lanenum', 'keepfastq', 'genome',
               'facility', 'machine', 'flowcell', 'rundate', 'url', 'runnumber')

  def __init__(self, libcode, lanenum, flowcell, rundate, runnumber=None,
               url=None, keepfastq=False, genome=None, 
               facility='EXT', machine='Unknown'):

    self.config    = Config()
    self.library   = Library.objects.get(code=libcode)
    self.lanenum   = lanenum
    self.keepfastq = keepfastq
    self.facility  = Facility.objects.get(code=facility)
    self.machine   = machine
    self.flowcell  = flowcell
    self.rundate   = rundate
    self.url       = url
    self.runnumber = runnumber

    if genome is None:
      self.genome = self.library.genome
    else:
      self.genome = Genome.objects.get(code=genome)

  def open_fastq(self, fastq):
    '''
    Return a file handle for a fastq file, transparently handling
    gzipped files.
    '''
    if is_zipped(fastq):
      from gzip import GzipFile
      gen = lambda x: GzipFile(x)
    else:
      gen = lambda x: open(x)
    return gen(fastq)

  def fastq_readlength(self, fastq):
    '''
    Guess the length of the reads in the fastq file. Assumes that the
    first read in the file is representative.
    '''
    # Currently just assumes that the second line is the first read, and
    # that it is representative.
    LOGGER.info("Finding read length from fastq file %s...", fastq)
    rlen = None
    with self.open_fastq(fastq) as reader:
      for _num in range(2):
        line = reader.next()
      rlen = len(line.rstrip('\n'))

    return rlen

  def fastq_readcount(self, fastq):
    '''
    Count the number of reads in the fastq file.
    '''
    LOGGER.info("Counting reads in fastq file %s...", fastq)
    lcount = 0
    with self.open_fastq(fastq) as reader:
      for line in reader:
        lcount += 1

    return lcount/4

  def create_unsaved_lane(self, fastqs):
    '''
    Gather statistics from the fastq files and generate a Lane object
    not yet saved to the database.
    '''
    status = Status.objects.get(code='complete')
    machine = Machine.objects.get(code__iexact=self.machine)

    # Don't save this to the db just yet.
    lane = Lane(library    = self.library,
                facility   = self.facility,
                lanenum    = self.lanenum,
                machine    = machine,
                flowcell   = self.flowcell,
                rundate    = self.rundate,
                summaryurl = self.url,
                runnumber  = self.runnumber,
                status     = status,
                flowlane   = 0)

    if len(fastqs) == 1:
      lane.paired = False
    elif len(fastqs) == 2:
      lane.paired = True
    else:
      raise ValueError("Can only process either one or two fastq files per lane.")

    # Assumes the first fastq is representative.
    lane.readlength   = self.fastq_readlength(fastqs[0])
    lane.reads        = self.fastq_readcount(fastqs[0])
    lane.passedpf     = lane.reads # This is a bit of an assumption FIXME?

    lanefiles = dict()
    if self.keepfastq:
      fqtype = Filetype.objects.get(code='fq')
      for fastq in fastqs:
        checksum = checksum_file(fastq)

        if is_zipped(fastq):
          fname = re.sub('%s$' % self.config.gzsuffix, '', fastq)
        else:
          fname = fastq
            
        try:
          lfile    = Lanefile.objects.get(lane=lane, filename=fname)
        except Lanefile.DoesNotExist:
          lfile    = Lanefile(lane     = lane,
                              filename = fname,
                              checksum = checksum,
                              filetype = fqtype)
          fastq = self._check_file_zipped(fastq, lfile)
          lanefiles[fastq] = lfile

    # Returns only the new lanefiles.
    return (lane, lanefiles)

  def create_unsaved_alignment(self, bam, lane):
    '''
    Gather statistics from the bam file, and generate a bed file and a
    new Alignment object unsaved to the database.
    '''
    bamtype = Filetype.objects.get(code='bam')
    bedtype = Filetype.objects.get(code='bed')

    # Generate a bed file from the bam.
    chrom_sizes = os.path.join(self.config.genomesizedir, self.genome.code + ".fa.length")
    if not os.path.exists(chrom_sizes):
      LOGGER.warning("Unable to find chromosome sizes file %s. BED file reads will be untrimmed.",
                     chrom_sizes)
      chrom_sizes = None
    bam2bed = BamToBedConverter(chrom_sizes=chrom_sizes)
    bambase = os.path.splitext(bam)[0]
    bed_fn  = bambase + bedtype.suffix
    beds = bam2bed.convert(bam, bed_fn)
    if len(beds) == 1:
      bed = beds[0]
    else:
      raise ValueError("Unexpected results from BAM to BED conversion")

    LOGGER.info("Counting reads in generated bed file %s...", bed)
    (mapped, unique) = count_reads(bed)

    # We don't save this yet because we're not currently within a
    # transaction.
    aln = Alignment(lane        = lane,
                    genome      = self.genome,
                    total_reads = lane.total_passedpf,
                    mapped      = mapped,
                    munique     = unique)

    # Create bam Alnfile.
    LOGGER.info("Checksumming bam file %s...", bam)
    checksum = checksum_file(bam)
    bamobj   = Alnfile(alignment=aln,
                       filetype=bamtype,
                       filename=bam, # casually assume no gzippng FIXME
                       checksum=checksum)

    # Create bed Alnfile
    LOGGER.info("Checksumming bed file %s...", bed)
    checksum = checksum_file(bed)
    bedobj   = Alnfile(alignment=aln,
                       filetype=bedtype,
                       filename=bed, # only just created, so not zipped
                       checksum=checksum)

    bamkey = self._check_file_zipped(bam, bamobj)
    bedkey = self._check_file_zipped(bed, bedobj)

    alnfiles = { bamkey : bamobj,
                 bedkey : bedobj }

    return (aln, alnfiles)

  def add(self, bam, fastqs=None, progname='bwa', progvers=None):
    '''
    Main entry point for the class.
    '''
    try:

      # A pre-existing lane is left almost unmolested.
      lane = Lane.objects.get(library  = self.library,
                              facility = self.facility,
                              lanenum  = self.lanenum)

    except Lane.DoesNotExist:

      # Creation of a new lane parses statistics from the fastq
      # file(s), runs fastqc and, if desired, will store the fastq file
      # itself in the repository.
      if fastqs is None:
        raise ValueError("Cannot create a new lane in the repository without"
                         + " the fastq files from which to harvest metadata.")
      (lane, lanefiles) = self.create_unsaved_lane(fastqs)
      self._save_lane_to_database(lane, lanefiles)

      if self.keepfastq:
        fastqs = None # Use the fastqs now stored in the repository.

      # Note: this code doesn't understand updating pre-existing reports.
      with LaneFastQCReport(target=lane, fastqs=fastqs, path=os.environ['PATH']) as fastqc:
        fastqc.insert_into_repository() # database transaction

    # Alignments are always appended to the lane; multiple alignments
    # may be added in this way.
    (aln, alnfiles) = self.create_unsaved_alignment(bam, lane)
    lane.mapped = aln.mapped
    lane.save()

    self._save_aln_to_database(aln, alnfiles, progname, progvers)

  def _check_file_zipped(self, fname, fobj):
    # Logging currently handled by the utilities module.
    zipped = is_zipped(fname)
    if fobj.filetype.gzip and not zipped:
      fname = rezip_file(fname, overwrite=True)
    elif not fobj.filetype.gzip and zipped:
      fname = unzip_file(fname, overwrite=True)
    return fname

  @transaction.atomic
  def _save_lane_to_database(self, lane, lanefiles):

    # Finally, save everything to database. We save before trying to
    # move the files because this step will raise an exception if
    # we're trying to overwrite a pre-existing file.
    lane.save()
    for lfile in lanefiles.values():
      lfile.lane = lane # ensure lane_id has been set.
      lfile.save()

    # Move fastq files into position, if that's a thing we're doing.
    if self.keepfastq:
      for (fastq, lfile) in lanefiles.iteritems():

        # Handle whether or not we're zipping fastq files.
        fastq = self._check_file_zipped(fastq, lfile) # defensive coding; this should already be set.
        dest  = lfile.repository_file_path
        destdir = os.path.dirname(dest)
        if not os.path.exists(destdir):
          os.makedirs(destdir)
        move(fastq, dest)
        set_file_permissions(self.config.group, dest)

  @transaction.atomic
  def _save_aln_to_database(self, aln, alnfiles, progname, progvers):
    # Handle the alignment.
    aln.save()
    for alf in alnfiles.values():
      alf.alignment = aln # ensure alignment_id has been set.
      alf.save()
      
    if progvers == None:
      alignerinfo = ProgramSummary(progname,
                                   ssh_host=self.config.cluster,
                                   ssh_user=self.config.clusteruser,
                                   ssh_path=self.config.clusterpath,
                                   ssh_port=self.config.clusterport)
      progname = alignerinfo.program
      progvers = alignerinfo.version

    try:
      program = Program.objects.get(program=progname,
                                    version=progvers,
                                    current=True)
    except Program.DoesNotExist, _err:
      raise StandardError("Unable to find current program in database: %s %s"
                          % (progname, progvers))

    DataProvenance.objects.create(program      = program,
                                  parameters   = '',
                                  rank_index   = 1,
                                  data_process = aln)

    for (fname, fobj) in alnfiles.iteritems():
      fname = self._check_file_zipped(fname, fobj) # defensive coding; this should already be set.
      dest  = fobj.repository_file_path
      destdir = os.path.dirname(dest)
      if not os.path.exists(destdir):
        os.makedirs(destdir)
      move(fname, dest)
      set_file_permissions(self.config.group, dest)

################################################################################
        
if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(\
    description='Add external data files to libraries already loaded in the repository.')

  PARSER.add_argument('-l', '--library', dest='library', type=str, required=True,
                      help='The library code with which the data should be associated.')

  PARSER.add_argument('-n', '--lanenum', dest='lanenum', type=int, default=1,
                      help='The number of the lane to store (or update) in the repository.'
                      + ' By default the script will add new alignments to lane number 1.'
                      + ' If this lane is already present in the database then the input'
                      + ' fastq file(s) will be ignored.')

  PARSER.add_argument('-b', '--bam', dest='bam', type=str, required=True,
                      help='The name of the bam file to load into the repository.')

  PARSER.add_argument('fastqs', metavar='<fastq files>', type=str, nargs='+',
                      help='The name of the fastq file or files to process.'
                      + ' The files may be a single or a pair of fastq'
                      + ' files. These files are not stored'
                      + ' in the repository by default (see --keep-fastq).')

  PARSER.add_argument('-k', '--keep-fastq', dest='keepfastq', action='store_true',
                      help='Whether or not to store the fastq file(s) in the'
                      + ' repository (default is to discard them).')
  
  PARSER.add_argument('-g', '--genome', dest='genome', type=str, required=False,
                      help='The genome code used for the alignment. The default is'
                      + ' to use the genome originally linked with the library.')

  PARSER.add_argument('--machine', dest='machine', type=str, default='Unknown',
                      help='The machine name of the sequencer.')

  PARSER.add_argument('--facility', dest='facility', type=str, default='EXT',
                      help='The facility code under which the lane should be categorised (default=EXT).')

  PARSER.add_argument('--flowcell', dest='flowcell', type=str, required=True,
                      help='The ID of the flowcell used in sequencing. This may be'
                      + ' set to the experiment accession if the flowcell ID is not available.')

  PARSER.add_argument('--rundate', dest='rundate', type=str, required=True,
                      help='The date of the sequencing run (may be the data submission date). YYYY-MM-DD.')

  PARSER.add_argument('--url', dest='url', type=str, required=False,
                      help='The experiment url in GEO/ArrayExpress or wherever.')

  PARSER.add_argument('--sample', dest='sample', type=str, required=False,
                      help='The accession code for the sample (e.g. GSM number in GEO).')

  PARSER.add_argument('--program', dest='program', type=str, default='bwa',
                      help='The name of the aligner used to generate the bam file. The'
                      + ' version will be that found in the pipeline cluster $PATH')

  PARSER.add_argument('--program-version', dest='progvers', type=str, default=None,
                      help='The version of the aligner program used. The default behaviour'
                      + ' is to connect to our cluster and find this out automatically;'
                      + ' this option is provided for those all-too-frequent occasions'
                      + ' nowadays when the cluster is unavailable.')

  ARGS = PARSER.parse_args()

  HND = ExternalDataHandler(libcode   = ARGS.library,
                            lanenum   = ARGS.lanenum,
                            facility  = ARGS.facility,
                            genome    = ARGS.genome,
                            keepfastq = ARGS.keepfastq,
                            machine   = ARGS.machine,
                            flowcell  = ARGS.flowcell,
                            rundate   = ARGS.rundate,
                            url       = ARGS.url,
                            runnumber = ARGS.sample)

  HND.add(bam      = ARGS.bam,
          fastqs   = ARGS.fastqs,
          progname = ARGS.program,
          progvers = ARGS.progvers)
