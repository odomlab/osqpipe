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

'''Script to update the adapters in the database given a tab-delimited
file with library codes in one column and adapter names in the
second.'''

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Library, Adapter

def update_adapters(filename):

  with open(filename) as infh:

    for line in infh:
      
      (lib, adpt) = line.split()

      libObj  = Library.objects.search_by_name(lib)
      adptObj = Adapter.objects.get(code=adpt)

      print("Updating library %s with adapter %s" % (lib, adpt))
      libObj.adapter = adptObj

      libObj.save()


if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Update adapters for libraries in the database.')

  PARSER.add_argument('-f', '--file', dest='file', type=str, required=True,
                      help='The input tab-delimited association file.')

  ARGS = PARSER.parse_args()

  update_adapters(ARGS.file)
