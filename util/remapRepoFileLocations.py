#!/usr/bin/env python

'''Quick script to automate the reorganisation of our repository filesystem.'''

# Currently just a testing skeleton, the commented lines below
# represent the actual working part.

import os
from shutil import move

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
