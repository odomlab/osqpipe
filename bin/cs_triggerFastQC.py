#!/usr/bin/env python

"""Creates FastQC repor for either a stand alone fastq file(s) or for a file from Odom repository for repository."""
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

from osqpipe.models import Lane, Lanefile, Alnfile
from osqpipe.pipeline.config import Config

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

def get_file_from_repository(conf, code, facility, lanenum, ftype):

  """Fetches files of a type from repository given code, facility and lane number. Returns list of files. NB! Note that function doed not currently work for QC files"""

  fns = []

  # loads file suffixes to file types information from Filetype table.
  try:
    lane = Lane.objects.get(library__code=code,
                            facility__code=facility,
                            lanenum=lanenum)
  except Lane.DoesNotExist, err:
    sys.exit("Failed to load lane '%s_%s%02d'.  Quitting." % (code, facility, lanenum))

  fns = Lanefile.objects.filter(lane=lane, filetype__code=ftype)
  fns = [ x for x in fns if os.path.exists(x.repository_file_path) ]
  
  # if file of type not found among lane files, search from alnfiles
  if not fns:

    fns = Alnfile.objects.filter(alignment__lane=lane,
                                 filetype__code=ftype)
    fns = [ x for x in fns if os.path.exists(x.repository_file_path) ]
    
  return fns

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


def unzipped_name(fn):
  """Returns fn without .gz suffix"""
  if os.path.splitext(fn)[1] == ".gz":
    fn = fn[0:-3]
  return fn

def compute_fast_qcforRepository(code, facility, replicate, dir, fastqc, java):
  
  """Computes FastQC report given the file in repository for repository, or if dir is defined not for repository but to dir"""  
  ftype = 'fq'

  conf = Config()

  fns = get_file_from_repository(conf, code, facility, replicate, ftype)
  if fns:
    ## Check if fastqc report has already been computed
    nofastqc = []
    for fn in fns:
      # get year from fastq file dir - not elegant! NOTE - NOW BROKEN by the repo reorg. FIXME.
      year = date.today().year
      matchObj = re.match('.*/(\d\d\d\d)/.*', fn, re.M|re.I)
      if matchObj:
        year = "%s" % matchObj.group(1)
      # check if report target dir exists
      if not dir:
        dir = os.path.join(conf.repositorydir, 'qcfile', "%s"%year, code)
      if not os.path.exists(dir):
        os.makedirs(dir)
      if not os.path.exists(dir):
        print "Failed to create FastQC files destination directory \"%s\". Quiting." % dir
        sys.exit(1)
      # check if FastQC report in target dir already exists
      unzippedfn = unzipped_name(fn)
      unzippedfn = os.path.basename(unzippedfn)
      qcfile = os.path.join(dir, unzippedfn + "_fastqc")
      if not os.path.exists(qcfile):
        nofastqc.append(fn)
      else:
        print "FastQC report in %s. Exiting." % (qcfile)
    ## Compute report if not disk.
    if nofastqc:
      compute_fast_qc(nofastqc, dir, fastqc, java)
    ## Insert computed fastqc files back to repository
  else:
    print "Failed to find files for Code=%s,  Facility=%s, Replicate=%s" % (code, facility, replicate)
    # sys.exit(1)

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
          compute_fast_qcforRepository(cols[0], cols[1], cols[2], "", fastqc, java)
          os.system("date")
        else:
          print "Failed to parse line \"%s\". Skipping." % line
      except ValueError:
        print "Failed to parse line \"%s\". Skipping." % line
  else:
    print "\nFile \"%s\" not accessible or missing! Quitting." % args.batch         
else:
  if args.code and args.facility and args.replicate:
    compute_fast_qcforRepository(args.code, args.facility, args.replicate, "", fastqc, java)
  else:
    print "\nUse either --input for a fastq file or specify values for all of the following: --code/--facility/--replicate; but not both.\n"
    parser.print_help()
    sys.exit(1)
