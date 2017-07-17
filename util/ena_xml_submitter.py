#!/usr/bin/env python
#
# $id$

import os
import sys
from osqpipe.pipeline.ena import enaXMLuploader
            
if __name__ == '__main__':
    
    from argparse import ArgumentParser
        
    PARSER = ArgumentParser(description='A script for uploading XMLs to ENA.')
    
    PARSER.add_argument('-t', '--xml_type', dest='xml_type', type=str, required=True, default=None,
                        help='XML type. E.g. SAMPLE or EXPERIMENT or RUN or ANALYSIS.')

    PARSER.add_argument('--validate', dest='validate_only', action='store_true', required=False, default=True,
                        help='Submit XML for validation only.')

    PARSER.add_argument('-c', '--ena_credenials_file', dest='credentials_file', type=str, required=True, default=None,
                        help='ENA credentials file is a tab delimited file where the first column is either \'username\' or \'password\' followed by values in the second.')
    
    PARSER_FILE_GROUP = PARSER.add_mutually_exclusive_group(required=True)
    
    PARSER_FILE_GROUP.add_argument('-d', '--submission_dir', dest='submission_dir', type=str, required=False, default=None,
                        help='Submission directory containing file md5s and XMLs for upload.')

    PARSER_FILE_GROUP.add_argument('-f', '--xml_file', dest='xml_file', type=str, required=False, default=None,
                        help='xml to be submitted. Note that the submitter expects also submission xml with same prefix in the file system. E.g. for 12345-sample.xml 12345-submission.xml is expected.')

    xml_types = ['SAMPLE','EXPERIMENT', 'RUN', 'ANALYSIS']
    
    ARGS = PARSER.parse_args()

    ena_loader = enaXMLuploader(ARGS.credentials_file)
    
    if ARGS.xml_type not in xml_types:
        print "xml type \"%s\" not recognized! Following types are supported: %s. Exiting!" % (ARGS.xml_type, ",".join(xml_types))
        sys.exit(1)

    if ARGS.xml_file is not None:
        if not os.path.isfile(ARGS.xml_file):
            print "%s not found! Exiting!" % ARGS.xml_file
        else:
            submission_xml = ARGS.xml_file[0:ARGS.xml_file.rindex('-')] + '-submission.xml'
            if os.path.isfile(submission_xml):
                accession = ena_loader.upload_xml(submission_xml,ARGS.xml_file,ARGS.xml_type, validate_only=ARGS.validate_only)
            else:
                print "%s not found! Exiting!" % submission_xml
    else:
        ena_loader.submit_files_in_directory(ARGS.submission_dir, ARGS.xml_type, ARGS.validate_only)
