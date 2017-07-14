#!/bin/env python

import os
import sys

from ExternalRecordManager import ExternalRecordManager

if __name__ == '__main__':
  
  from argparse import ArgumentParser
  PARSER = ArgumentParser(description='Check and add external records to lanes, libraries and samples.')

  PARSER_OPTIONS = PARSER.add_mutually_exclusive_group(required=True)  
  PARSER_OPTIONS.add_argument('--ena_receipt_file', dest='ena_receipt', type=str,
                      help='Name of the ENA submission receipt XML file.')
  PARSER_OPTIONS.add_argument('--library','-l', dest='library', type=str, help='Library code.')
  PARSER_OPTIONS.add_argument('--sample','-s', dest='sample_name', type=str, help='Sample name.')
  PARSER_OPTIONS.add_argument('--lane_id', dest='lane_id', type=int, help='Lane ID.')
  PARSER_OPTIONS.add_argument('--project', dest='project_name', type=str, help='Project name. Adds external record to all lanes part of the project.')  
  PARSER.add_argument('--accession','-a', dest='accession', type=str, help='External accession.', required=False, default=None)
  
  PARSER.add_argument('--public', dest='public', action='store_true', help='Data has been released to public.', required=False, default=False)
  PARSER.add_argument('--release_date', dest='release_date', type=str, help='Date of data release (YYYY-MM-DD).', required=False)
  PARSER.add_argument('--repository_name', dest='repository_name', type=str, help='Repository name. E.g. \'EBI ENA\', \'EBI BioStudies\'.', required=True)
  PARSER.add_argument('-t', dest='test_only', action='store_true',
                      help='Test if external record exists.')
  
  ARGS = PARSER.parse_args()
  
  if ARGS.ena_receipt:
    erm = ExternalRecordManager()
    erm.parse_xml_receipt(ARGS.ena_receipt)
    erm.check_external_record_in_obj()
    if not ARGS.test_only:
      erm.add_external_record_to_obj()
    sys.exit(0)

  if ARGS.accession is None:
    sys.stdout.stderr("Accession required!")
    sys.exit(1)

  if ARGS.project_name:
    # find lanes part of the project:
    try:
      lanes = Lane.objects.filter(library__projects__name=ARGS.project_name)
    except Lane.DoesNotExist:
      sys.stdout.write("No lanes associated with project \'%s\'.", ARGS.project_name)
      sys.exit(1)
    for lane in lanes:
      erm = ExternalRecordManager(entry_type='lane', accession = ARGS.accession, release_date = ARGS.release_date, is_public = ARGS.public, repository_name=ARGS.repository_name)
      erm.add_lane_obj(lane.id)
      erm.check_external_record_in_obj()
      if not ARGS.test_only:
        erm.add_external_record_to_obj()
    sys.exit(0)
  
  if ARGS.library:
    if ARGS.sample_name:
      erm = ExternalRecordManager(entry_type='sample', accession = ARGS.accession, release_date = ARGS.release_date, is_public = ARGS.public, repository_name=ARGS.repository_name)
      erm.add_sample_obj_library(ARGS.library)
    else:
      erm = ExternalRecordManager(entry_type='library', accession = ARGS.accession, release_date = ARGS.release_date, is_public = ARGS.public, repository_name=ARGS.repository_name)
      erm.add_library_obj(ARGS.library)
  elif ARGS.sample_name and not ARGS.library:
    erm = ExternalRecordManager(entry_type='sample', accession = ARGS.accession, release_date = ARGS.release_date, is_public = ARGS.public, repository_name=ARGS.repository_name)
    erm.add_sample_obj(ARGS.sample_name)
  elif ARGS.lane_id:
    erm = ExternalRecordManager(entry_type='lane', accession = ARGS.accession, release_date = ARGS.release_date, is_public = ARGS.public, repository_name=ARGS.repository_name)
    erm.add_lane_obj(ARGS.lane_id)

  erm.check_external_record_in_obj()
  if not ARGS.test_only:
    erm.add_external_record_to_obj()
    sys.exit(0)
