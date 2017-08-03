#!/usr/bin/env python

'''
Very simple script to check that all the raw data from a flowcell
has been downloaded and loaded into the repository, and that each Lane
object is linked to at least one alignment. This script will also flag
up cases where alignment rates are lower than expected.
'''

from logging import INFO

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane, Library
from osqpipe.pipeline.flowcell import FlowCellQuery
from osqpipe.pipeline.fastq_aligner import FastqBwaAligner, FastqTophatAligner

from osqutil.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)

###############################################################################

def confirm_complete(flowcell, bwa_algorithm,
                     resubmit_alignments=False):
  '''
  Given a flowcell ID, confirm completion of the pipeline. Returns
  True for complete, False for incomplete.
  '''
  # Suppress some of the more verbose output here, as it's likely to
  # no longer be of interest.
  fcq = FlowCellQuery(flowcell_id = flowcell, quiet = True)

  complete = True
  for ( flowlane, libcodeset ) in fcq.lane_library.iteritems():
    for libcode in sorted(list(libcodeset)):
      try:
        lane = Lane.objects.get(library__code=libcode,
                                flowcell=flowcell,
                                flowlane=flowlane)
      except Lane.DoesNotExist:
        LOGGER.warning("Lane not found in repository: %s %s:%s",
                       libcode, flowcell, flowlane)
        complete = False
        continue
      if lane.lanefile_set.filter(filetype__code__in=('fq','tar')).count() == 0:
        LOGGER.warning("Lane has no fastq files or tarballs: %s %s:%s",
                       libcode, flowcell, flowlane)
        complete = False
        continue
      if lane.alignment_set.count() == 0:
        LOGGER.warning("Lane has no alignments: %s %s:%s",
                       libcode, flowcell, flowlane)
        complete = False
        if resubmit_alignments:
          library = Library.objects.get(code=libcode)
          if library.libtype.code == 'rnaseq':
            aligner = FastqTophatAligner(samplename=library.sample.name)
          else:
            aligner = FastqBwaAligner(samplename=library.sample.name,
                                      bwa_algorithm=bwa_algorithm)

          aligner.align(library  = libcode,
                        facility = lane.facility.code,
                        lanenum  = lane.lanenum,
                        genome   = library.genome.code)

        continue
      for aln in lane.alignment_set.all():
        if aln.alnfile_set.filter(filetype__code='bam').count() == 0:
          LOGGER.warning("Alignment has no bam files: %s %s:%s (%s)",
                         libcode, flowcell, flowlane, aln.genome)
          complete = False
        if aln.mapped_percent < 50:
          LOGGER.warning("Alignment has low mapping rate:"
                         + " %s %s:%s (%s): %d%%",
                         libcode, flowcell, flowlane, aln.genome,
                         aln.mapped_percent)

  return complete

###############################################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(\
    description='Confirm completion of FlowCell processing.')

  PARSER.add_argument('flowCell', metavar='<flowCellID>', type=str,
                      help='The flow cell ID (required).')

  PARSER.add_argument('-r', '--resubmit', dest='resubmit', action='store_true',
                      help='Indicates that missing alignments should be'
                      + ' resubmitted against the standard configured genome'
                      + ' for the library. Output bam files will be written'
                      + ' to the current working directory.')

  PARSER.add_argument('--algorithm', type=str, dest='algorithm',
                      choices=('aln', 'mem'),
                      help='If --resubmit is specified, the bwa algorithm to'
                      + ' use (aln or mem). The default behaviour is to pick'
                      + ' the algorithm based on the read length in the fastq'
                      + ' files.')

  ARGS = PARSER.parse_args()

  if confirm_complete(flowcell            = ARGS.flowCell,
                      resubmit_alignments = ARGS.resubmit,
                      bwa_algorithm       = ARGS.algorithm):
    LOGGER.info("Flowcell processing complete.")
  else:
    LOGGER.warning("Flowcell processing not yet complete.")
