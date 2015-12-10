#!/usr/bin/env python

'''Very simple script to find files in the directories under the
current working directory, search the repository for the corresponding
Alnfile, and update the date column if found. Usage: change directory
to e.g. repository/bam, and just run this script.'''

import os
import time
import re

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Alnfile

gz_re = re.compile('\.gz$')

for root, dirs, files in os.walk('.', followlinks=True):

  for filename in files:

    dbfile = gz_re.sub('', filename)

    try:
      fobj = Alnfile.objects.get(filename=dbfile)
    except Alnfile.DoesNotExist, err:
      continue

    path = os.path.join(root, filename)

    st = os.stat(path)

    modtime = time.gmtime(st.st_mtime)

    date = "%04d-%02d-%02d" % (modtime.tm_year, modtime.tm_mon, modtime.tm_mday)

    fobj.date = date

    print "updating file %s with date %s" % (filename, date)

    fobj.save()
