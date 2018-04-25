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

"""
Creates FastQC report for either a stand alone fastq file(s) or for a
file from Odom repository for repository.
"""
__author__ = "Margus Lukk"
__date__ = "09 Nov 2012"
__version__ = "0.1"
__credits__ = ""

#
# Dependencies:
# - fastqc
# - java
# - half dependency on contraction of repository qc file directories.
#
# Todo:
# - in case of repository related FastQC report, check (from DB) if report already exists
# - in case of repository related FastQC report, add report files to repository
#


import sys
import os
import os.path
import re
import argparse
from datetime import date

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane, Lanefile, Alnfile
from osqutil.config import Config
from osqpipe.pipeline.laneqc import LaneFastQCReport

def compute_fast_qc(fns, target, fastqc, java, threads=2):

  """Executes fastqc report generation"""

  # consider having other options in here
  if len(fns) < threads:
    threads = len(fns)
  cmd = "%s -q -t %s -o %s -j %s" % (fastqc, threads, target, java)
  for fn in fns:
    cmd = cmd + " %s" % (fn)
  print "%s" % cmd
  os.system(cmd)

def check_files_or_dirs(list, message, type="file"):

  """Checks if files or dirs in the list exist. Returns first existing element in a list and a log message."""
  
  pos = 0
  log = ""
  for e in list:
    if e:
      if type == "dir" and not os.path.exists(e):
        os.makedirs(e)
      if os.path.exists(e):
        return (e, log)
      elif pos < len(list)-1:
        log = log + "%s not accessible or found! Trying alternative.\n" % e
      else: 
        log = log + "%s not accessible or found! %s\n" % (e, message)
  return ("", log)

def compute_fast_qcforRepository(code, facility, replicate):
  
  """
  Computes and stores a FastQC report for a lane in the
  repository. This function will raise an exception if the lane
  already has a fastqc report. Note that this code links into the
  standard pipeline report generator and so will correctly produce
  PDFs as well as the regular report files.
  """
  conf = Config()

  lane = Lane.objects.get(library__code=code, facility__code=facility, lanenum=replicate)

  if lane.laneqc_set.filter(provenance__program__program='fastqc').count() == 0:
    with LaneFastQCReport(target=lane, path=conf.hostpath) as qcrep:
      qcrep.insert_into_repository()
  else:
    raise StandardError("Lane already has a FastQC report.")

################## M A I N ########################

# Set some default variables
JAVA='/home/fnc-odompipe/software/external/bin/java'
FASTQC='/home/fnc-odompipe/software/external/bin/fastqc'
dir = "."

parser = argparse.ArgumentParser(
  description='Executes FastQC report generation.')
parser.add_argument('-q', '--fastqc', dest='fastqc', required=False, type=str, help='Location of fastqc.')
parser.add_argument('-j', '--java', dest='java', required=False, type=str, help='Location of java.')
parser.add_argument('-t', '--targetdir', dest='targetdir', required=False, type=str, help='Target directory for FastQC files.')
parser.add_argument('-c', '--code', dest='code', required=False, type=str, help='The library code. E.g. do123.')
parser.add_argument('-f', '--facility', dest='facility', required=False, type=str, help='The name of the facility. E.g. \"CRI\".')
parser.add_argument('-r', '--replicate', dest='replicate', required=False, type=str, help='The replicate number. E.g. 1.')
parser.add_argument('-b', '--batch', dest='batch', required=False, type=str, help='Name of a file containing tab separated code, facility and replicate on each row. Use the option alone.')
parser.add_argument('file', nargs='*', type=str, help='Fastq file name. Can be compressed by gzip.')
args = parser.parse_args()

# Check if java and fastq installations exist
(java, jm) = check_files_or_dirs([args.java, JAVA], "Specify location of java with --java.")
if jm:
  print "%s" % jm
if not java:
  parser.print_help()
  sys.exit(1)
(fastqc, fm) = check_files_or_dirs([args.fastqc, FASTQC], "Specify location of fastq with --fastqc.")
if fm:
  print "%s" % fm
if not fastqc:
  parser.print_help()
  sys.exit(1)

# Decide if to run with --input file option or with --code/--facility/--replicate options either in batch or entered with flags
if args.file:
  if args.code or args.facility or args.replicate:
    print "\nUse either --input for a fastq file or specify values for all of the following: --code/--facility/--replicate; but not both.\n"
    parser.print_help()
    sys.exit(1)
  else:
    (dir, dm) = check_files_or_dirs([args.targetdir], "Use --targetdir.")
    if dm:
      print "%s" %dm
    if not dir:
      parser.print_help()
      sys.exit(1)
    compute_fast_qc(args.file, dir, fastqc, java)
# batch file with Code\tFacility\tReplicate or single Code\tFacility\tReplicate entry?
elif args.batch:
  if os.path.exists(args.batch):
    for line in open(args.batch):
      line = line.rstrip()
      try:
        cols = line.split("\t")
        if cols[0] and cols[1] and cols[2]:
          os.system("date")
          print "Processing \"%s\" \"%s\" \"%s\"..." % (cols[0], cols[1], cols[2])
          compute_fast_qcforRepository(cols[0], cols[1], cols[2])
          os.system("date")
        else:
          print "Failed to parse line \"%s\". Skipping." % line
      except ValueError:
        print "Failed to parse line \"%s\". Skipping." % line
  else:
    print "\nFile \"%s\" not accessible or missing! Quitting." % args.batch         
else:
  if args.code and args.facility and args.replicate:
    compute_fast_qcforRepository(args.code, args.facility, args.replicate)
  else:
    print "\nUse either --input for a fastq file or specify values for all of the following: --code/--facility/--replicate; but not both.\n"
    parser.print_help()
    sys.exit(1)
