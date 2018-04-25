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

'''Quick script to automate the reorganisation of our repository filesystem.'''

# Currently just a testing skeleton, the commented lines below
# represent the actual working part.

import os
from shutil import move

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Alnfile, Lanefile

# Keeping these as stub methods so that remap*FileLocations.py scripts
# don't cause errors with pylint.
def old_build_alnfile_path(alnfile):
  raise NotImplementedError()

def old_build_lanefile_path(lanefile):
  raise NotImplementedError()

def rename_files(files, fromFun):
  for fobj in files:

    old = fromFun(fobj)
    new = fobj.repository_file_path

    if old != new:
      if os.path.exists(old):
        d = os.path.dirname(new)
        if not os.path.exists(d):
          os.makedirs(d)
          print "Created directory %s" % (d,)
        print "Moving %s to %s" % (old, new)
        move(old, new)
      else:
#        print "MISSING FILE: %s" % (old,)
        raise OSError("MISSING FILE (gzip error?): %s" % (old,))


files = Alnfile.objects.all()
rename_files(files, old_build_alnfile_path)
    
files = Lanefile.objects.all()
rename_files(files, old_build_lanefile_path)
