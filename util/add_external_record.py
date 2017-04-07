#!/bin/env python

import os
import sys
from datetime import datetime
from subprocess import Popen, PIPE
import xml.etree.ElementTree as ET

import django
from osqpipe.models import Lane, Library, Sample, ExternalRecord, ExternalRepository
from osqpipe.models import Machine

from logging import INFO, DEBUG
from osqutil.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)

django.setup()

class ExternalRecordManager(object):
  
  def __init__(self, entry_type=None, accession = None, release_date = None, is_public = False, repository_name='EBI ENA'):

    self.receipt_type = None
    self.accession = None
    self.release_date = None
    self.obj_type = None
    self.obj = None
    self.repository_name = repository_name
    self.is_public = is_public
    
    if entry_type is not None:
      self.entry_type = entry_type
    if accession is not None:
      self.accession = accession
    if release_date is not None:
      self.release_date = release_date      

  def add_external_record_to_obj(self):
    '''Inserts external_record and links it with the object obj which is either lane, sample or library.'''
    
    release_date = datetime.strptime(self.release_date , '%Y-%m-%d')
    # 1. We may want to make release date optional in which case current date would be used.
    # 2. We may want to override is_public in case release_date is in past
    try:
      er = ExternalRecord.objects.get(accession__iexact=self.accession)
    except ExternalRecord.DoesNotExist:
      erep = ExternalRepository.objects.get(name__iexact=self.repository_name)
      er = ExternalRecord(accession=self.accession, repository=erep, is_public=self.is_public, release_date=self.release_date)
      er.save()

    self.obj.external_records.add(er)
    self.obj.save()
    
    print "%s linked to %s %s" % (self.accession, self.obj_type, self.obj)

  def check_external_record_in_obj(self):
    '''Check if object (lane, library or sample) has external record for accession in self.accession'''
    ers = None
    if self.obj_type == 'lane':
      try:
        ers = ExternalRecord.objects.filter(lanes=self.obj, repository__name='EBI ENA')
      except ExternalRecord.DoesNotExist:
        print "No external recors associated with lane %s." % self.obj
    if self.obj_type == 'library':
      try:
        ers = ExternalRecord.objects.filter(libraries=self.obj, repository__name='EBI ENA')
      except ExternalRecord.DoesNotExist:
        print "No external recors associated with library %s." % self.obj
    if self.obj_type == 'sample':
      try:
        ers = ExternalRecord.objects.filter(samples=self.obj, repository__name='EBI ENA')
      except ExternalRecord.DoesNotExist:
        print "No external recors associated with sample %s." % self.obj

    found = False
    if len(ers) == 0:
      print "No ENA entries found!"
      return 0
    else:
      for er in ers:
        if er.accession == self.accession:
          found = True
          print "%s accession %s already associated with %s %s." % (er.repository.name, er.accession, self.obj_type, self.obj)
          return 1

    if not found:
      print "No entry for accession %s found!" % self.accession
      return 0


  def add_lane_obj(self, lane_id):
    try:
      self.obj = Lane.objects.get(id=lane_id)
    except Lane.DoesNotExist:
      sys.stdout.write("No lane with id %d.", lane_id)
      sys.exit(1)
    self.obj_type = 'lane'

  def add_library_obj(self, code):
    try:
      self.obj = Library.objects.get(code=code)
    except Library.DoesNotExist:
      print "No library with code \"%s\"." % code
      sys.exit(1)
    self.obj_type = 'library'

  def add_sample_obj(self, sample_name):
    try:
      self.obj = Sample.objects.get(name=sample_name)
    except Sample.DoesNotExist:
      print "No sample with name \"%s\"." % sample_name
      sys.exit(1)
    self.obj_type = 'sample'


  def parse_xml_receipt(self, fname):
    '''Parse ENA submission receipt XML'''

    alias = None
    receipt_types = ['RUN','EXPERIMENT','SAMPLE']
    
    tree = ET.parse(fname)
    root = tree.getroot()
    # RECEIPT receiptDate="2017-01-24T11:14:48.479Z
    if root.tag == 'RECEIPT':
      self.release_date = root.get('receiptDate')[0:10]
    for child in root:
      if child.tag in receipt_types:
        self.receipt_type = child.tag
        self.accession = child.get('accession')
        alias = child.get('alias')
        # In case external accession has been specified, use that instead.
        # <EXT_ID accession="SAMEA53947168" type="biosample"/>
        for gchild in child:
          if gchild.tag == "EXT_ID":
            self.accession = gchild.get('accession')

    # Double check that we have all the values.
    failed = False
    if self.release_date is None:
      print "Failed to parse %s. Release_date is missing!" % (fname)
      failed = True
    if self.accession is None:
      print "Failed to parse %s. Acccession is missing!" % (fname)
      failed = True
    if self.receipt_type is None:
      print "Failed to parse %s. Receipt is of unknown type. (%s expected)." % (fname, " or ".join(receipt_types))
      failed = True
    if alias is None:
      print "Failed to parse %s. Failing to parse \"%s\" into anything meaningful." % (fname, alias)
      failed = True
    if failed:
      sys.exit(1)

    # Parse alias into obj and obj_type
    if self.receipt_type == 'RUN':
      # e.g. do9147_CRI_3_Flowcell=C9FMUANXX_Lane=5
      (code, facility, lanenum, flowcell_tmp, flowlane_tmp) = alias.split('_')
      (label, flowcell) = flowcell_tmp.split('=')
      (label, flowlane) = flowlane_tmp.split('=')
      try:
        self.obj = Lane.objects.get(library__code=code,facility__code=facility, lanenum=lanenum, flowcell=flowcell, flowlane=flowlane)
      except Lane.DoesNotExist:
        sys.stdout.write("No lane found for code=%s, facility=%s, lanenum=%d, flowcell=%s, flowlane=%s.", code, facility, lanenum, flowcell, flowlane)
        sys.exit(1)
      self.obj_type = 'lane'

    if self.receipt_type == 'EXPERIMENT':
      self.add_library_obj(code=alias)

    if self.receipt_type == 'SAMPLE':
      self.add_sample_obj(sample_name=alias)

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
  PARSER.add_argument('--accession','-a', dest='accession', type=str, help='External accession.', required=True)
  
  PARSER.add_argument('--public', dest='public', action='store_true', help='Data has been released to public.', required=False, default=False)
  PARSER.add_argument('--release_date', dest='release_date', type=str, help='Date of data release (YYYY-MM-DD).', required=False)
  PARSER.add_argument('--repository_name', dest='repository_name', type=str, help='Repository name. E.g. \'EBI ENA\', \'EBI BioStudies\'.', required=True)
  PARSER.add_argument('-t', dest='test_only', action='store_true',
                      help='Test only if external record has already been added.')
  
  ARGS = PARSER.parse_args()
  
  if ARGS.ena_receipt:
    erm = ExternalRecordManager()
    erm.parse_xml_receipt(ARGS.ena_receipt)
    sys.exit(0)

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
    erm = ExternalRecordManager(entry_type='library', accession = ARGS.accession, release_date = ARGS.release_date, is_public = ARGS.public, repository_name=ARGS.repository_name)
    erm.add_lane_obj(ARGS.library)

  if ARGS.sample_name:
    erm = ExternalRecordManager(entry_type='sample', accession = ARGS.accession, release_date = ARGS.release_date, is_public = ARGS.public, repository_name=ARGS.repository_name)
    erm.add_sample_obj(ARGS.sample_name)

  if ARGS.lane_id:
    erm = ExternalRecordManager(entry_type='lane', accession = ARGS.accession, release_date = ARGS.release_date, is_public = ARGS.public, repository_name=ARGS.repository_name)
    erm.add_lane_obj(ARGS.lane_id)

  erm.check_external_record_in_obj()
  if not ARGS.test_only:
    erm.add_external_record_to_obj()
