#!/usr/bin/env python
#
# $Id$

'''Generate a histogram summarising duplicate reads.'''

import sys
import os
import logging

from fastq import FastqIO
from histogram import Histogram

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()

LOG_LEVEL = logging.DEBUG
ISATTY = False
ADAPTER = 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCT'

################################################################################

def init():
  '''Set some global output options.'''
  global ISATTY
  if (os.isatty(sys.stderr.fileno())):
    ISATTY = True
  LOGGER.setLevel(LOG_LEVEL)
  logfmt = "[%(asctime)s] %(levelname)s : %(message)s"
  fmt = logging.Formatter(logfmt)
  hdlrErr = logging.StreamHandler(sys.stderr)
  hdlrErr.setFormatter(fmt)
  LOGGER.addHandler(hdlrErr)

def _init_reads():
  '''Initialise some variables.'''
  LOGGER.debug("loading reads...")
  reads = {}
  total = 0
  adap = 0
  return (reads, total, adap)

def _add_read(data, reads, total, adap):
  '''Compare a sequence against our nascent dict structure and add a
  count in the appropriate place. Also increment total and adapter
  counts as necessary.'''
  if data[0:len(ADAPTER)] == ADAPTER:
    adap += 1
  else:
    if data in reads:
      reads[data] += 1
    else:
      reads[data] = 1
    total += 1
    if total % 100000 == 0 and ISATTY:
      sys.stderr.write("%9d\r" % (total,))
  return (total, adap)

def load_reads_fq(fname):
  '''Read the fastq file into a dict keyed by sequence. Also return a
  count of adapter reads and a total read count.'''
  (reads, total, adap) = _init_reads()
  fdesc = FastqIO(fname)
  for seq in fdesc:
    data = seq.seq
    (total, adap) = _add_read(data, reads, total, adap)
  fdesc.close()
  return (reads, total, adap)

def make_histo(reads):
  '''Given a dict of reads keyed by sequence, create a Histogram
  object summarising the distribution of duplicated reads.'''
  LOGGER.debug("making histogram...")
  total = 0
  histo = Histogram()
  for read in reads.iterkeys():
    histo.add(reads[read])
    if reads[read] >= 1000:
      sys.stdout.write("%s\t%d\n" % (read, reads[read]))
    total += 1
    if total % 100000 == 0 and ISATTY:
      sys.stderr.write("%9d\r" % (total,))
  return histo

def count_dups(histo):
  '''Count the total duplicates from the output of make_histo.'''
  dups = 0
  for (num, count) in histo.iterBuckets():
    dups += (num - 1) * count
  return dups

def perc(num, denom):
  '''Simple percent calculator function.'''
  pcnt = int(((float(num) / float(denom)) * 1000)) / 10.0
  return "%4.1f%%" % (pcnt,)

def summarize(fname):
  '''Generate a histogram summarising duplicate reads.'''
  (reads, total, adap) = load_reads_fq(fname)
  histo = make_histo(reads)
  duplicates = count_dups(histo)
  sys.stdout.write("total %d\nduplicated %d (%s)\nadapter %d (%s)\n" % (
      total, duplicates, perc(duplicates, total),
      adap, perc(adap, total+adap)))
  histo.dump(sys.stdout, includeEmpty=False)  

################################################################################

if __name__ == '__main__':

  init()
  summarize(sys.argv[1])
