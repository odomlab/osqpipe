#!/usr/bin/env python
#
# $Id$

'''Utility functions used to generate quality plots.'''

import sys
import os
import time

from Gnuplot import Gnuplot

from pipeline.config import Config
CONFIG = Config()

############################################################
def plot_qv(base, means, stdevs):
  '''
  Generate a quality plot and return the URL to the output file.
  '''
  filename = "%s/%s" % (CONFIG.httptmpdir, base)
  url      = "%s/%s" % (CONFIG.httptmpurl, base)

  # There's almost certainly a prettier way to do this in the Gnuplot
  # package; this is a straight port from piping commands to the
  # gnuplot program:
  g = Gnuplot()
  g('set terminal png size 400,300')
  g("set output \"%s\"" % filename)
  g("plot [0:%d] [-1:41] \"-\" with errorbars notitle," % (len(means)+1,)
    + " \"-\" with linespoints notitle")
  for i in range(len(means)):
    g("%d %f %f" % (i+1, means[i], stdevs[i]))
  g("e")
  for i in range(len(means)):
    g("%d %f %f" % (i+1, means[i], stdevs[i]))
  g("e")
  g("set output")
  return url

def plot_pfqual_values(lane):
  '''Plot quality values for reads passing vendor check ("passed filter").'''
  base = "qualplot%d_%dpf.png" % (os.getpid(), int(time.time()))
  return plot_qv(base, lane.qualmeanpf, lane.qualstdevpf)

def plot_all_qual_values(lane):
  '''Plot quality values for all reads.'''
  base = "qualplot%d_%dnpf.png" % (os.getpid(), int(time.time()))
  return plot_qv(base, lane.qualmean, lane.qualstdev)

############################################################

