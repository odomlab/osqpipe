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

import os
import sys
import xml.etree.ElementTree as ET
from logging import INFO, DEBUG
from osqutil.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)


from ena_submit_library import enaXMLuploader

class EnaCancelXml(object):

    def __init__(self, accession=None, xml_template_dir=None, xml_target_dir=None):

        self.xml_target_dir = xml_target_dir
        if xml_target_dir is None:
            self.xml_target_dir = os.getcwd()

        self.submission_cancel_template = os.path.join(xml_template_dir, 'submission_cancel.xml')
        if not os.path.isfile(self.submission_cancel_template):
            LOGGER.error("XML template for submission cancels not found! File expected: %s" % self.submission_cancel_template)
            sys.exit(1)

    def write_cancel_xml(self, accession):
        '''Takes ENA accession and writes cancel XML for this accession.'''

        # Set full path for the output XML file
        xml_file = os.path.join(self.xml_target_dir, accession + "-submission_cancel.xml")        
        if os.path.isfile(xml_file):
            LOGGER.info("Submission cancel XML for \'%s\' already exists (%s)!", accession, xml_file)
            sys.exit(0)

        # Parse XML template and
        xml = ET.ElementTree()        
        xml.parse(self.submission_cancel_template)
        # Find CANCEL tag and modify target attribute to the accession
        submission = xml.find('SUBMISSION')
        submission.set('center_name', 'CANCER RESEARCH UK CAMBRIDGE INSTITUTE')
        actions = submission.find('ACTIONS')
        action = actions.find('ACTION')        
        cancel = action.find('CANCEL')
        cancel.set('target', accession)
        
        # Write submission XML
        xml.write(xml_file)

        return xml_file

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='ENA submission cancel utility.')
  PARSER.add_argument('-a', '--accession', dest='accession', type=str, required=True,
                                      help='Accession to be cancelled.')

  PARSER.add_argument('-x', '--xml_template_dir', dest='xml_template_dir', type=str, required=True,
                      help='Directory containing ENA short read data submission templates.')
  PARSER.add_argument('--xml_target_dir', dest='xml_target_dir', type=str, required=False, default="",
                      help='Directory where XML files should be created.')

  PARSER.add_argument('-v', '--validate_xml', dest='validate_xml', action='store_true', default=False,
                      help='Create XML files and validate against ENA test submission system but do not submit the data.')
  PARSER.add_argument('--submit_xml', dest='submit_xml', action='store_true', default=False,
                      help='Submit XMls to ENA.')

  PARSER.add_argument('-c', '--ena_credentials_file', dest='ena_credentials_file', type=str, required=True,
                      help='File containing ENA credentials with username on first line and password on second.')

  ARGS = PARSER.parse_args()
  # install submitter                                                                                                                                       
  ena_cancel = EnaCancelXml(xml_template_dir=ARGS.xml_template_dir, xml_target_dir=ARGS.xml_target_dir)
  submission_xml_filename = ena_cancel.write_cancel_xml(ARGS.accession)

  ena_loader = enaXMLuploader(ARGS.ena_credentials_file)
  if ARGS.validate_xml:
      ena_loader.upload_xml(submission_xml_filename, validate=True)
  if ARGS.submit_xml:
      ena_loader.upload_xml(submission_xml_filename, validate=False)
