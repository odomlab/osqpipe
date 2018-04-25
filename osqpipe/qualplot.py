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

'''Utility functions used to generate quality plots.'''

import sys
import os
import time

from Gnuplot import Gnuplot

from osqutil.config import Config
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
  # gnuplot program (Update FIXME seems like gnuplotlib is the way to
  # go, or move to something completely unrelated to gnuplot like
  # matplotlib or seaborn):
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

