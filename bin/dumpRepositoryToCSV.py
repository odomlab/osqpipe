#!/usr/bin/env python

'''
Script to automatically dump a core set of library and lane annotation
out to CSV format for sharing with collaborators (via e.g. Dropbox).
'''

from osqpipe.models import Lane

from osqpipe.pipeline.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

class RepositoryDumper(object):
  '''
  Class which dumps a specific set of metadata from all libraries and
  lanes in the repository to an output file.
  '''
  __slots__ = ('mapping', 'separator')

  def __init__(self, separator="\t"):
    # This mapping controls the output; the list contains 2-element
    # tuples which are (header_string, lambda function to retrieve row
    # contents).
    self.mapping = [ ('libcode',  lambda x: x.library.code),
                     ('tissue',   lambda x: x.library.tissue.name),
                     ('facility', lambda x: x.facility.code) ]
    self.separator = separator

  def header_string(self):
    '''
    Returns a string to be used as file header, ready for writing to
    output file.
    '''
    return self.separator.join([ elem[0] for elem in self.mapping ])

  def laneobj_to_string(self, lane):
    '''
    Returns a string representing the passed lane object, ready for
    writing to output file.
    '''
    return self.separator.join([ elem[1](lane) for elem in self.mapping ])

  def dump_to_file(self, outfile):
    '''
    Dumps everything to the specified output file.
    '''
    with open(outfile, 'w') as outfh:
      outfh.write(self.header_string())
      lanes = Lane.objects.all()
      for lane in lanes:
        lanestr = self.laneobj_to_string(lane)
        outfh.write(lanestr)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='Dump library metadata to CSV output.')

  PARSER.add_argument('-o', '--output', dest='output', type=str, required=True,
                      help='The name of the output file.')

  ARGS = PARSER.parse_args()

  DUMPER = RepositoryDumper()

  DUMPER.dump_to_file(ARGS.output)
