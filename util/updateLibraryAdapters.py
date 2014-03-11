#!/usr/bin/env python

'''Script to update the adapters in the database given a tab-delimited
file with library codes in one column and adapter names in the
second.'''

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
