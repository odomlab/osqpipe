#!/bin/env python
#
# $id$

import os
import sys

from osqpipe.pipeline.ena import EnaSubmitter

if __name__ == '__main__':

  from argparse import ArgumentParser

  # Required arguments
  
  PARSER = ArgumentParser(description='Data submission to pre-exisgting ENA study.')

  # Define what needs to be submitted:
  # - whole project (based on project name),
  # - individual lane (based on lane id),
  # - a set of lanes listed in a file (lane ids),
  # - or a library (library code)
  PARSER_P_OR_LANE_GROUP = PARSER.add_mutually_exclusive_group(required=True)  
  PARSER_P_OR_LANE_GROUP.add_argument('-p', '--project', dest='project', type=str,
                                      help='The name of the project in chipseq repository.')
  PARSER_P_OR_LANE_GROUP.add_argument('-l', '--lane_id', dest='lane_id', type=str,
                                      help='Lane ID in chipseq repository.')
  PARSER_P_OR_LANE_GROUP.add_argument('-L', '--lane_ids_file', dest='lane_id_file', type=str,
                                      help='A file containing lane IDs for which the data should be submitted in the first column.')
  PARSER_P_OR_LANE_GROUP.add_argument('--libcode', dest='libcode', type=str,
                                      help='Library code.')

  # Define study at ENA where the data will be added by providing one of the following:
  # - providing submission XML of the study
  # - study alias
  # - study accessionOB
  PARSER_STUDY_GROUP = PARSER.add_mutually_exclusive_group(required=True)
  PARSER_STUDY_GROUP.add_argument('-s', '--study_xml', dest='study_xml', type=str,
                      help='Name of the ENA study XML file. Needed for extraction of ENA study alias. Data will be added to study with this alias at ENA.')
  PARSER_STUDY_GROUP.add_argument('-S', '--study_alias', dest='study_alias', type=str,
                                  help='Study alias as submitted to ENA. Data will be added to study with this alias at ENA.')
  PARSER_STUDY_GROUP.add_argument('-a', '--study_accession', dest='study_accession', type=str,
                                  help='Study accession in ENA. Data will be added to study with this accession at ENA.')   
  # Define submission level:
  # - upload all files (and create md5 sums if not present in the file path)
  # - create XMLs
  # - validate XMLs
  # - submit XMLs
  PARSER.add_argument('--upload_files', dest='upload_files', action='store_true', default=False,
                      help='Uploads short read files related to submission')
  PARSER.add_argument('-v', '--validate_xml', dest='validate_xml', action='store_true', default=False,
                      help='Create XML files and validate against ENA test submission system but do not submit the data.')
  PARSER.add_argument('--submit_xml', dest='submit_xml', action='store_true', default=False,
                      help='Submit XMls to ENA.')

  PARSER.add_argument('-x', '--xml_template_dir', dest='xml_template_dir', type=str, required=True,
                      help='Directory containing ENA short read data submission templates.')

  PARSER.add_argument('-c', '--ena_credentials_file', dest='ena_credentials_file', type=str, required=True,
                      help='File containing ENA credentials with username on first line and password on second.')
  # Optional arguments:
  # - a local directory for storing XMLs created in submission process
  # - release date
  # - validate XMLs by submitting them to ENA validation server.
  PARSER.add_argument('--xml_target_dir', dest='target_dir', type=str, required=False, default="",
                      help='Directory where XML files should be created.')
  
  PARSER.add_argument('-d', '--release_date', dest='release_date', type=str, required=False, default=None,
                      help='Release date of the lane/experiment. Defaults to a date 2 years from now.')
  
  ARGS = PARSER.parse_args()
  # install submitter
  ena = EnaSubmitter(study_alias=ARGS.study_alias, study_xml=ARGS.study_xml, xml_template_dir=ARGS.xml_template_dir, credentials_file=ARGS.ena_credentials_file, release_date_str=ARGS.release_date, validate_xml=ARGS.validate_xml, upload_files=ARGS.upload_files, submit_xml=ARGS.submit_xml, target_dir=ARGS.target_dir)
  
  if ARGS.project:
    ena.submit_project(project=ARGS.project, )

  if ARGS.lane_id:
    # e.g. lane_id=75696
    ena.submit_lane_with_id(ARGS.lane_id)
  
  if ARGS.lane_id_file:
    ena.submit_lanes_from_file(fn=ARGS.lane_id_file)

  if ARGS.libcode:
    ena.submit_library(libcode=ARGS.libcode)
