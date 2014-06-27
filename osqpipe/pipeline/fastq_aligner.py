#!/usr/bin/env python
#
# $Id$

'''Code to rerun the alignment of a given library's fastq file against
a genome registered in the repository.'''

import os
import logging

from ..models import Library, Genome, Filetype
from config import Config
from bwa_runner import BwaClusterJobSubmitter, BwaDesktopJobSubmitter

from setup_logs import configure_logging
LOGGER = configure_logging()

class FastqAligner(object):

  '''Abstract base class used to handle the repository queries prior
  to running any alignment.'''

  __slots__ = ('conf', 'test_mode', 'finaldir')

  def __init__(self, test_mode, finaldir=None):
    self.conf = Config()
    self.test_mode = test_mode
    if test_mode:
      LOGGER.setLevel(logging.DEBUG)
    else:
      LOGGER.setLevel(logging.INFO)

    # Default to saving the output in the current working directory.
    if finaldir is None:
      finaldir = os.path.realpath(os.getcwd())
    self.finaldir = finaldir

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

  def align(self, library, facility, lanenum, genome,
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

    # We specifically want fastq ("fq") files.
    wanted = Filetype.objects.get(code='fq')

    LOGGER.debug("Retrieving fastq files for library from database.")
    lanestoprocess = lanes.all()
    if lanenum is not None:
      lanestoprocess = lanestoprocess.filter(lanenum=lanenum)
    if facility is not None:
      lanestoprocess = lanestoprocess.filter(facility__code=facility)

    for lane in lanestoprocess:
      files = []
      lanecount = 0
      for fobj in [ x for x in lane.lanefile_set.all()
                    if x.filetype == wanted ]:
        files += [ fobj.repository_file_path ]
        lanecount += 1
      if lanecount > 2:
        raise ValueError("More than one fastq file found in database for lane.")

      if len(files) == 0:
        LOGGER.warning(
          "No fastq files found in database for lane %s (%s). Skipping.",
          lane, library)
        return
      elif len(files) > 2:  # up to 2 (for paired-end reads).
        raise ValueError("More than two fastq files found in database"
                         + " for lane %d (%s)." % (lane.lane, library))

      self._call_aligner(files, gobj, destnames=destnames,
                        nocc=nocc, cleanup=(not nocleanup))

  def align_standalone(self, filepaths, genome, destnames=None,
                       nocc=None, nocleanup=False):

    '''Align fastq file(s) against a genome recorded in the
    repository. Does not require the fastq files themselves to be
    present in the repository.'''

    gobj = self._retrieve_genome(genome)
    self._call_aligner(filepaths, gobj, destnames=destnames,
                      nocc=nocc, cleanup=(not nocleanup))

class FastqBwaAligner(FastqAligner):

  '''Class used to handle alignments using the BWA external
  binary. This is currently coded to use BWA running on a cluster via
  the BwaClusterJobSubmitter class.'''

  def _call_aligner(self, filepaths, genome, destnames=None,
                    nocc=None, cleanup=True, *args, **kwargs):

    '''Method used to dispatch cs_runBwaWithSplit.py processes on the
    cluster.'''

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
                    finaldir=self.finaldir,
                    num_threads=num_threads)
    bsub.submit(filenames=filepaths, auto_requeue=False,
                destnames=destnames,
                is_paired=paired, cleanup=cleanup, nocc=nocc)

    LOGGER.info("Jobs submitted.")

