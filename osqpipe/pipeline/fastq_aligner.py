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

'''Code to rerun the alignment of a given library's fastq file against
a genome registered in the repository.'''

import os

from ..models import Library, Genome, Filetype
from osqutil.config import Config
from osqutil.utilities import determine_readlength
from .bwa_runner import BwaClusterJobSubmitter, BwaDesktopJobSubmitter, TophatClusterJobSubmitter, StarClusterJobSubmitter

from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging()

class FastqAligner(object):

  '''Abstract base class used to handle the repository queries prior
  to running any alignment.'''

  __slots__ = ('conf', 'test_mode', 'finaldir', 'samplename')

  def __init__(self, test_mode=False, finaldir=None, samplename=None):
    self.conf = Config()
    self.test_mode = test_mode
    if test_mode:
      LOGGER.setLevel(DEBUG)
    else:
      LOGGER.setLevel(INFO)

    # Default to saving the output in the current working directory.
    if finaldir is None:
      finaldir = os.path.realpath(os.getcwd())
    self.finaldir = finaldir
    self.samplename = samplename

  def _call_aligner(self, filepaths, genome, *args, **kwargs):

    '''Stub method used in subclasses to call out to external
    alignment program.'''

    raise NotImplementedError()

  @staticmethod
  def _retrieve_genome(genome):
    '''
    Very simple wrapper method to retrieve a genome from the database.
    '''
    LOGGER.debug("Querying database for genome.")
    try:
      gobj = Genome.objects.get(code=genome)
    except Genome.DoesNotExist, _err:
      raise StandardError("Genome code not found in database: %s" % (genome,))
    return gobj

  def align(self, library, genome, facility=None, lanenum=None,
            destnames=None, nocc=None, nocleanup=False):

    '''Core method to retrieve all fastq files for a given library,
    confirm that the requested genome is in the database, and then
    dispatch a suitable aligner on the retrieved files.'''

    gobj = self._retrieve_genome(genome)

    LOGGER.debug("Querying database for library.")
    try:
      lib = Library.objects.search_by_name(library)
    except Library.DoesNotExist, _err:
      raise StandardError("Library code not found in database: %s" % (library,))

    # If PolIII/TFIIIC library, set number of non-unique reads as in
    # self.conf.nonuniquereads
    if lib.factor and lib.factor.name in self.conf.reallocation_factors:
      nocc = self.conf.nonuniquereads

    LOGGER.debug("Retrieving sequencing lanes for library from database.")
    lanes = lib.lane_set
    if lanes.count() == 0:
      raise StandardError("No sequencing lanes found for library: %s"
                          % (library,))

    if lib.libtype.code in ('smrnaseq', 'ripsmrnaseq'):
      # SmallRNA-seq uses fasta files from the reaper/tally pipeline.
      wanted_filetype = 'fa'
    else:
      # For every thing else, we want fastq ("fq") files.
      wanted_filetype = 'fq'

    LOGGER.debug("Retrieving %s files for library from database.", wanted_filetype)
    lanestoprocess = lanes.all()
    if lanenum is not None:
      lanestoprocess = lanestoprocess.filter(lanenum=lanenum)
    if facility is not None:
      lanestoprocess = lanestoprocess.filter(facility__code=facility)

    for lane in lanestoprocess:
      files = []
      for fobj in lane.lanefile_set.filter(filetype__code__iexact=wanted_filetype):
        files += [ fobj.repository_file_path ]

      if len(files) == 0:
        LOGGER.warning(
          "No %s files found in database for lane %s. Skipping.",
          wanted_filetype, lane)
        continue

      if len(files) > 1 and not lane.paired:
        raise ValueError("Multiple %s files found for single-ended lane %s."
                         % (wanted_filetype, lane))

      elif len(files) > 2:  # (paired-end reads only).
        raise ValueError(("More than two %s files found in database"
                         + " for lane %s.") % (wanted_filetype, lane))

      self._call_aligner(files, gobj, destnames=destnames,
                        nocc=nocc, cleanup=(not nocleanup))

  def align_standalone(self, filepaths, genome, destnames=None,
                       nocc=None, nocleanup=False, nosplit=False, rcp=None, lcp=None, fileshost=None):

    '''Align fastq file(s) against a genome recorded in the
    repository. Does not require the fastq files themselves to be
    present in the repository.'''

    gobj = self._retrieve_genome(genome)
    self._call_aligner(filepaths, gobj, destnames=destnames,
                       nocc=nocc, cleanup=(not nocleanup), nosplit=nosplit, rcp=rcp, lcp=lcp, fileshost=fileshost)

class FastqBwaAligner(FastqAligner):
  '''
  Class used to handle alignments using the BWA external
  binary. This is currently coded to use BWA running on a cluster via
  the BwaClusterJobSubmitter class.
  '''
  def __init__(self, bwa_algorithm=None, *args, **kwargs):
    assert(bwa_algorithm in (None, 'aln', 'mem'))
    super(FastqBwaAligner, self).__init__(*args, **kwargs)
    self.bwa_algorithm = bwa_algorithm

  def _choose_bwa_algorithm(self, filepaths):
    '''
    Method determines read length and sets the bwa algorithm to 'mem'
    if 70bp or greater, else 'aln'.
    '''
    # If we want to modify the default algorithm selection, this is a
    # good place to do so.
    rlen = determine_readlength(filepaths[0])
    LOGGER.debug("FASTQ read length: %d", rlen)
    return 'mem' if rlen >= 70 else 'aln'

  def _call_aligner(self, filepaths, genome, destnames=None,
                    nocc=None, cleanup=True, nosplit=False, rcp=None, lcp=None, fileshost=None, *args, **kwargs):
    '''
    Method used to dispatch cs_runBwaWithSplit.py processes on the
    cluster.
    '''
    if len(filepaths) == 1:
      LOGGER.info("Launching single-end sequencing alignment.")
      paired = False
    elif len(filepaths) == 2:
      LOGGER.info("Launching paired-end sequencing alignment.")
      paired = True
    else:
      raise ValueError("Wrong number of files passed to aligner.")

    if destnames is None:
      destnames = [os.path.basename(x) for x in filepaths]

    if self.bwa_algorithm is None:
      self.bwa_algorithm = self._choose_bwa_algorithm(filepaths)
    LOGGER.info("BWA algorithm chosen: %s", self.bwa_algorithm)

    # We use the alternative alignment host mechanism if it's been
    # configured; otherwise we default to using the standard
    # cluster-based mechanism.
    num_threads = 1
    try:
      althost = self.conf.althost
      assert(althost != '')
    except AttributeError, _err:
      althost = None
    if althost is None:
      jobclass = BwaClusterJobSubmitter
    else:
      jobclass = BwaDesktopJobSubmitter
      num_threads = self.conf.althostnumthreads

    # Build path to genome index on the alignment host/cluster.
    genome_path = jobclass.build_genome_index_path(genome)

    bsub = jobclass(test_mode=self.test_mode,
                    genome=genome_path,
                    samplename=self.samplename,
                    finaldir=self.finaldir,
                    num_threads=num_threads)
    bsub.submit(filenames=filepaths, auto_requeue=False,
                destnames=destnames, bwa_algorithm=self.bwa_algorithm,
                is_paired=paired, cleanup=cleanup, nocc=nocc, nosplit=nosplit, rcp=rcp, lcp=lcp, fileshost=fileshost)

    LOGGER.info("Jobs submitted.")

class FastqTophatAligner(FastqAligner):
  '''
  Class used to handle alignments using the tophat2 external
  binary. This is currently coded to use tophat2 running on a cluster via
  the TophatClusterJobSubmitter class.
  '''
  def _call_aligner(self, filepaths, genome, destnames=None,
                    nocc=None, cleanup=True, *args, **kwargs):
    '''
    Method used to dispatch cs_runTophatWithSplit.py processes on the
    cluster.
    '''
    if nocc is not None:
      LOGGER.warning("Unsupported nocc argument passed to FastqTophatAligner.")

    if len(filepaths) == 1:
      LOGGER.info("Launching single-end sequencing alignment.")
      paired = False
    elif len(filepaths) == 2:
      LOGGER.info("Launching paired-end sequencing alignment.")
      paired = True
    else:
      raise ValueError("Wrong number of files passed to aligner.")

    if destnames is None:
      destnames = [os.path.basename(x) for x in filepaths]

    # The alternative alignment host mechanism is currently
    # unsupported for tophat2 alignments (for lack of development time).
    num_threads = 1
    try:
      althost = self.conf.althost
      assert(althost != '')
    except AttributeError, _err:
      althost = None
    if althost is None:
      jobclass = TophatClusterJobSubmitter
    else:
      raise ValueError("Currently unable to run tophat2 jobs via the"
                       + " alternative alignment host mechanism.")

    # Build path to genome index on the alignment host/cluster.
    genome_path = jobclass.build_genome_index_path(genome)

    bsub = jobclass(test_mode=self.test_mode,
                    genome=genome_path,
                    samplename=self.samplename,
                    finaldir=self.finaldir,
                    num_threads=num_threads)
    bsub.submit(filenames=filepaths, auto_requeue=False,
                destnames=destnames,
                is_paired=paired, cleanup=cleanup)

    LOGGER.info("Jobs submitted.")

class FastqStarAligner(FastqAligner):
  '''
  Class used to handle alignments using the STAR external
  binary. This is currently coded to use STAR running on a cluster via
  the StarClusterJobSubmitter class.
  '''
  def _call_aligner(self, filepaths, genome, destnames=None,
                    nocc=None, cleanup=True, *args, **kwargs):
    '''
    Method used to dispatch cs_runStarWithSplit.py processes on the
    cluster.
    '''
    if nocc is not None:
      LOGGER.warning("Unsupported nocc argument passed to FastqStarAligner.")

    if len(filepaths) == 1:
      LOGGER.info("Launching single-end sequencing alignment.")
      paired = False
    elif len(filepaths) == 2:
      LOGGER.info("Launching paired-end sequencing alignment.")
      paired = True
    else:
      raise ValueError("Wrong number of files passed to aligner.")

    if destnames is None:
      destnames = [os.path.basename(x) for x in filepaths]

    # The alternative alignment host mechanism is currently
    # unsupported for STAR alignments (for lack of development time).
    num_threads = 1
    try:
      althost = self.conf.althost
      assert(althost != '')
    except AttributeError, _err:
      althost = None
    if althost is None:
      jobclass = StarClusterJobSubmitter
    else:
      raise ValueError("Currently unable to run STAR jobs via the"
                       + " alternative alignment host mechanism.")

    # Build path to genome index on the alignment host/cluster.
    genome_path = jobclass.build_genome_index_path(genome)

    bsub = jobclass(test_mode=self.test_mode,
                    genome=genome_path,
                    samplename=self.samplename,
                    finaldir=self.finaldir,
                    num_threads=num_threads, aligner='star')
    bsub.submit(filenames=filepaths, auto_requeue=False,
                destnames=destnames,
                is_paired=paired, cleanup=cleanup)

    LOGGER.info("Jobs submitted.")
