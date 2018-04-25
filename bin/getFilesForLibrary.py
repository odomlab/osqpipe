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

'''List files available for a given sample. The script is currently
configured to search only for FASTQ files but this is easily
changed.'''

import re
from urlparse import urlparse

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqutil.config import Config
from osqpipe.pipeline.upstream_lims import runs_containing_samples,\
    get_lims_run_details

ARTIFACT_LABEL = {
  'FASTQ':     re.compile(r'Read \d FASTQ$'),
  'FASTQC':    re.compile(r'FASTQC Lane Report$'),
  'MGA':       re.compile(r'MGA Lane Report$'),
  'FASTQ_MD5': re.compile(r'FASTQ MD5 Checksums'),
  }

###############################################################################

def get_files(libcode, filetype='FASTQ'):
  '''List sequencing data files available for a given LIMS sample ID.'''
  conf = Config()

  ## If we know how to get the file type, use the standard artifact
  ## name; otherwise just treat filetype as a regular expression.
  label = ARTIFACT_LABEL.setdefault(filetype, re.compile(filetype))

  path_re = re.compile(r'(.*)/([^\/]+)$')

  root = runs_containing_samples(conf.lims_rest_uri, libcode)

  count = 0
  for run_elem in root.findall('./run'):
    run_id = run_elem.find('./runFolder').text

    print "Run ID: %s\n" % run_id
    
    run_root = get_lims_run_details(conf.lims_rest_uri, run_id)

    for lib_elem in run_root.findall("./run/flowcell/library/sample[name='%s']/.." % libcode):
      for file_elem in lib_elem.findall('./file'):
        name = file_elem.find('./artifactName').text
        if label.search(name):
          url = urlparse(file_elem.find('./url').text)
          pathbits = path_re.match(url.path)
          if pathbits:
            print ("host: %s:%d\npath: %s\nfile: %s\n"
                   % (url.hostname, url.port, pathbits.group(1), pathbits.group(2)))
            count += 1
          else:
            raise ValueError("Unexpected URL path structure: %s" % url.path)  

  if count == 0:
    print "Nothing found."

###############################################################################

if __name__ == '__main__':
  from argparse import ArgumentParser
  PARSER = ArgumentParser(description="Generate a listing of LIMS files available for a given library.")

  PARSER.add_argument('libcode', metavar='<libcode>', type=str,
                      help='The library code name (required).')

  PARSER.add_argument('-t', '--type', dest='filetype', type=str, default='FASTQ',
                      help='Specify the file type to use in the query (default: FASTQ).')

  ARGS = PARSER.parse_args()

  get_files(ARGS.libcode, ARGS.filetype)
