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
from pipes import quote
from subprocess import CalledProcessError

from osqutil.utilities import call_subprocess, bash_quote
from osqutil.config import Config

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Alnfile, Lanefile
from remapRepoFileLocations import old_build_alnfile_path, old_build_lanefile_path

SERVERHOST='10.20.13.6'
SERVERUSER='admin'
SERVERDIR='/share/MD0_DATA/backups/uk-cri-lsrv01'
IDKEY='/home/fnc-odompipe/.ssh/id_dsa'

CONFIG = Config()

def ssh_command(command):

  cmd = " ".join(('ssh', '-i', IDKEY,
                  "%s@%s" % (SERVERUSER, SERVERHOST)) + command)
  return(cmd)

def rename_files(files, fromFun):
  for fobj in files:

    old = fromFun(fobj)
    old = "/".join((SERVERDIR, old))
    new = fobj.repository_file_path
    new = "/".join((SERVERDIR, new))

    if old != new:
      d = os.path.dirname(new)
        
#      print "Creating directory %s" % (d,)
      cmd = ssh_command(('mkdir', '-p', quote(d)))
      call_subprocess(cmd, shell=True, path=CONFIG.hostpath)

#      print "Moving %s to %s" % (old, new)
      old = bash_quote(old)
      new = bash_quote(new)
      cmd = ssh_command(('mv', quote(old), quote(new)))
      try:
        call_subprocess(cmd, shell=True, path=CONFIG.hostpath)
      except CalledProcessError, err:
        print "Warning: move failed for file %s: %s" % (old, err)

files = Alnfile.objects.all()
rename_files(files, old_build_alnfile_path)
    
files = Lanefile.objects.all()
rename_files(files, old_build_lanefile_path)
