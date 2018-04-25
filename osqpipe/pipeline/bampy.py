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
Code tools used to manipulate bam files from within python. This is
intended to be a pure-python implementation for improved
maintainability.
'''

import os
import re
import pysam
from contextlib import contextmanager
from logging import INFO
from osqutil.setup_logs import configure_logging
LOGGER = configure_logging('bampy', level=INFO)

class Bamfile(pysam.AlignmentFile):

  '''
  Subclass of pysam.AlignmentFile with a few convenience methods of our own.
  '''

  # Note that while it might be nice to have a Bamfile class which
  # correctly auto-generates its own index, this is not going to be
  # easily maintainable due to the implementation of
  # pysam.AlignmentFile. Rather than coding ourselves into knots, we have
  # opted for a context manager approach (see open_bamfile).

  def read_refname(self, read):
    '''
    Return the name of the reference sequence to which a given read
    has been aligned.
    '''
    return self.getrname(read.tid)

  def basename(self):
    '''
    Return the name of the underlying filename, minus its extension.
    '''
    return os.path.splitext(self.filename)[0]

  def has_index_file(self):
    '''
    Simply establish whether a .bai index file is present for this bam
    file.
    '''
    return os.path.exists("%s.bai" % self.filename)

  def split_by_chromosome(self, pattern=None):
    '''
    Split reads by chromosome, writing out new bam files to disk and
    returning new Bamfile objects. Takes an optional pattern argument
    to enable filtering, e.g. for those genomes with many small
    scaffolds.
    '''
    bamlist = []
    if pattern is not None:
        pattern = re.compile(pattern)
    for refname in self.references:
      if pattern is not None:
        if not pattern.match(refname):
          continue
      LOGGER.info("Writing out reads for reference %s...", str(refname))
      self.reset()
      outfile = "%s_%s.bam" % (self.basename(), str(refname))
      out = pysam.AlignmentFile(outfile, 'wb', header=self.header)
      for read in self.fetch(reference=refname):
        out.write(read)
      out.reset()
      bamlist.append(Bamfile(filename=outfile))
      
    return bamlist

  def expand_tallied_reads(self, outfile, keep_multi=False):
    '''
    Expand a bam file containing reaper/tally output, using the format
    from the Odom Lab sequencing pipeline, into a bam file in which
    all reads are represented by their own individual records. If the
    keep_multi flag is set as True, multi-mapping reads will be kept.
    Note that this can be misleading as multi-mapping reads are no
    longer distributed randomly between mapped loci.
    '''
    LOGGER.info("Expanding tallied bam file %s", self.filename)
    self.reset()
    out = pysam.AlignmentFile(outfile, 'wb', header=self.header)
    for read in self.fetch(until_eof=True): # All reads in file, in order.

      # If mapping_quality is zero, we assume this means it's
      # multi-mapping (this is the bwa default).
      if not keep_multi and read.mapping_quality == 0:
        continue
      readparts = read.query_name.split(':')
      if len(readparts) == 2:
        tally = readparts[1].split('_')
        if len(tally) == 2 and tally[0] == 'count':
          count = int(tally[1])
        else:
          raise StandardError(\
            "Found read ID not matching expected tally output format: %s" % read.query_name)
      else:
        raise StandardError(\
          "Found read ID not matching expected tally output format: %s" % read.query_name)
      base_name = readparts[0] + ':read_'
      for n in range(count):
        read.query_name = base_name + str(n+1)
        out.write(read)

    out.close()
    return Bamfile(filename=outfile)

@contextmanager
def open_bamfile(filename, index=True, *args, **kwargs):
  '''
  Context manager function which automatically creates .bai index
  files if they're missing.
  '''
  if not os.path.exists("%s.bai" % filename) and index:
    LOGGER.info("Indexing bam file %s...", filename)
    pysam.index(filename)
  bam = Bamfile(filename=filename, *args, **kwargs)
  try:
    yield bam
  finally:
    bam.close()

if __name__ == '__main__':

  # This is purely here as test code; do not use this as a production script!
  import sys
  with open_bamfile(filename=sys.argv[1]) as bam:
    bam.split_by_chromosome()
  
