#!/usr/bin/env python
#
# $Id$

""" Fetch MGA reports from Solexa LIMS """
__author__ = "Margus Lukk"
__date__ = "22 Oct 2012"
__version__ = "0.1"
__credits__ = "Exploits cs_* libraries implemented by Gordon Brown."

## Dependencies:
#
# - Password free login to target scp server
# - scp
# - sed
# - /home/fnc-odompipe/software/external/bin/xsltproc
# - /home/fnc-odompipe/software/external/bin/wkhtmltopdf-amd64

import sys
import getopt

from osqpipe.pipeline.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

from osqpipe.pipeline.fetch_mga import fetch_mga

TEST_MODE = False

def usage(fname):
  """Prints usage """
  print ""
  print "Program: %s (fetches MGA reports from Solexa LIMS)" % (fname,)
  print "Version: %s" % __version__
  print "Usage:   %s [-t/--test] <flowcell> <lane> <dest> [filename prefix] " % (fname,)
  print "Options: -t/--test          run script in test mode"
  print "         <flowcell>         flowcell ID"
  print "         <lane>           lane number"
  print "         <dest>             location for saving the files"
  print "         [filename prefix]  prefix for local file names"
  print "Example: %s D0VJGACXX 5 /home/username/ do3333_liver_CRI01.mga" % (fname,)
  print ""

def process_cmdline(argv):
  """Process command line arguments"""
  global TEST_MODE
  try:
    (opts, args) = getopt.gnu_getopt(argv[1:], 't', ['test'])
  except getopt.GetoptError:
    usage(sys.argv[0])
    sys.exit("Incorrect script arguments.")
  for (key, _val) in opts:
    if key in ('-t', '--test'):
      TEST_MODE = True
  if len(args) < 3 or len(args) > 4:
    usage(sys.argv[0])
    sys.exit("Incorrect script arguments.")
  if len(args) == 3:
    nameprefix = "%s_%s_mga" % (args[0], args[1])
    args.append(nameprefix)
  return args

################ M A I N ######################################

if __name__ == '__main__':

  (FLOWCELL, FLOWLANE, DESTINATION, NAMEPREFIX) = process_cmdline(sys.argv)
  fetch_mga(FLOWCELL, FLOWLANE, DESTINATION,
           NAMEPREFIX)
