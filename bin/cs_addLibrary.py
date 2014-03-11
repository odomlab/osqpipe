#!/usr/bin/env python
#
# $Id$

'''Given a set of metadata on the command-line, create a new row in
the library repository table.'''

from osqpipe.pipeline.library import LibraryHandler

from osqpipe.pipeline.config import Config

from osqpipe.pipeline.setup_logs import configure_logging
LOGGER = configure_logging()
CONFIG = Config()

######################################################################

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Add library metadata to the repository.')

  PARSER.add_argument('project', type=str, default=CONFIG.defaultproject,
                      help='The project name.')

  PARSER.add_argument('libtype', type=str,
                      help='The type of library (e.g. chipseq, rnaseq).')

  PARSER.add_argument('label', type=str,
                      help='The library code (do#).')

  PARSER.add_argument('genome', type=str,
                      help='The genome to align the sequence reads against.')

  PARSER.add_argument('tissue', type=str,
                      help='The tissue from which the library was derived.')

  PARSER.add_argument('--strain', dest='strain', type=str,
                      help='The strain from which the library was derived.')

  PARSER.add_argument('--description', dest='description', type=str,
                      help='The individual from which the library was derived.')

  PARSER.add_argument('--factor', dest='factor', type=str,
                      help='The (typically ChIPseq) factor'
                      + ' being assayed by the library.')

  PARSER.add_argument('--antibody', dest='antibody', type=str,
                      help='The actual antibody used in preparing'
                      + ' the library (i.e. ChIPseq).')

  PARSER.add_argument('--lot_number', dest='lot_number', type=str,
                      help='The lot number for the antibody.')

  PARSER.add_argument('--individual', dest='individual', type=str,
                      help='The individual from which the library was derived.')

  PARSER.add_argument('--chipsample', dest='chipsample', type=str,
                      help='The experiment code associated with this library.')

  PARSER.add_argument('--barcode', dest='barcode', type=str,
                      help='The barcode associated with this library.')

  PARSER.add_argument('--linkerset', dest='linkerset', type=str,
                      help='The linker set used in library construction.')

  PARSER.add_argument('--paired', dest='paired', type=str, default='n',
                      help='Whether the library is for paired-end sequencing'
                      + ' or not. This should be a boolean value as recognised'
                      + ' by postgresql (e.g., yes/no).')

  PARSER.add_argument('--adapter', dest='adapter', type=str,
                      help='The adapter used in library construction.')

  ARGS = PARSER.parse_args()

  OPTS = vars(ARGS)

  REQD = ('project', 'libtype', 'label', 'genome', 'tissue')
  OPTS = dict( (x, y) for (x, y) in OPTS.iteritems()
               if x not in REQD and y is not None )

  HND = LibraryHandler()

  HND.add(projcodes = [ARGS.project],
          libtype = ARGS.libtype,
          code    = ARGS.label,
          genome  = ARGS.genome,
          tissue  = ARGS.tissue,
          opts    = OPTS)

