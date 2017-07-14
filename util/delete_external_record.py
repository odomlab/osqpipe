#!/bin/env python
import os
import sys

from ExternalRecordManager import ExternalRecordManager

if __name__ == '__main__':
  
  from argparse import ArgumentParser
  PARSER = ArgumentParser(description='Check and add external records to lanes, libraries and samples.')

  PARSER.add_argument('--ena_receipt_file', dest='ena_receipt', type=str,help='Name of the ENA submission receipt XML file.')
  # PARSER.add_argument('--accession','-a', dest='accession', type=str, help='External accession.', required=False, default=None)
  PARSER.add_argument('--libcode', dest='libcode', type=str, help='Library code.', required=False)
  PARSER.add_argument('--sample_name', dest='sample_name', type=str, help='Sample name.', required=False)

  PARSER.add_argument('-t', dest='test_only', action='store_true', default=False,
                      help='Test only if external record has already been added.')
  
  ARGS = PARSER.parse_args()
  
  if ARGS.ena_receipt:
    erm = ExternalRecordManager(code=ARGS.libcode)
    erm.remove_external_obj_from_receipt(ARGS.ena_receipt, test_only=ARGS.test_only)
