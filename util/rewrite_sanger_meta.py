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
Trivial but probably useful script, which takes a sequencing metadata
file from the Sanger pipeline and rewrites all sample IDs as defined
by a second, mapping file. The mapping file is assumed to consist of
two columns of tab-delimited values; column 1 contains the IDs to be
changed, and column 2 contains the replacement IDs.
'''

# N.B. SCRIPT IN DEVELOPMENT. CURRENTLY AWAITING DATA FROM BIANCA FOR
# TESTING.

import os
import re

def read_mapping(mapfile):
  '''
  Read the mapping file into a dict and return it.
  '''
  mapping = {}
  with open(mapfile) as infh:
    for line in infh:
      (orig, new) = [x.strip() for x in line.split("\t")]
      mapping[orig] = new

  return mapping

def rewrite(metafile, mapfile):
  '''
  Rewrite the metadata file. Renames the original file to
  <original>.old. Note that running the script twice will therefore
  potentially destroy the original data.
  '''
  mapping = read_mapping(mapfile)
  attr_re = re.compile(r"^attribute:\s*(.*?)\s*$")
  newfile = "%s.rewritten" % metafile
  with open(newfile, 'wb') as outfh:
    with open(metafile, 'rb') as infh:
      line = infh.readline()
      while line != '':
        outline = line
        attr_match = attr_re.match(line)
        if attr_match:
          if attr_match.group(1) in ('sample', 'library'):
            outfh.write(outline)
            outline = infh.readline()
            for (orig, new) in mapping.iteritems():
              padded_orig = r' %s\b' % orig
              if re.search(padded_orig, outline):
                outline = re.sub(padded_orig, ' %s' % new, outline)
        outfh.write(outline)
        line = infh.readline()
  os.rename(metafile, "%s.old" % metafile)
  os.rename(newfile, metafile)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(
    description='Script to remap sample IDs in Sanger metadata files.')

  PARSER.add_argument('-f', '--file', dest='file', type=str, required=True,
                      help='The metadata file to rewrite')

  PARSER.add_argument('-m', '--mapping', dest='mapping', type=str,
                      required=True, help='The mapping file to use.')

  ARGS = PARSER.parse_args()

  rewrite(ARGS.file, ARGS.mapping)
