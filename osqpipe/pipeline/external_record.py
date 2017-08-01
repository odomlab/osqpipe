#!/bin/env python
#
# $Id$

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
    '''A class for managing insertions (and in future if needed also deletions) of External Records'''
    def __init__(self, entry_type=None, accession = None, release_date = None, is_public = False, repository_name=None):
      
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

        if self.release_date is None:
            today = datetime.date.today()
            release_date = today.replace(year=today.year + 2)
        else:
            release_date = datetime.datetime.strptime(self.release_date , '%Y-%m-%d')
        # 1. We may want to make release date optional in which case current date would be used.
        # 2. We may want to override is_public in case release_date is in past
        try:
            er = ExternalRecord.objects.get(accession__iexact=self.accession)
        except ExternalRecord.DoesNotExist:
            erep = ExternalRepository.objects.get(name__iexact=self.repository_name)
            er = ExternalRecord(accession=self.accession, repository=erep, is_public=self.is_public, release_date=release_date)
            er.save()

        self.obj.external_records.add(er)
        self.obj.save()

        LOGGER.info("%s linked to %s %s", self.accession, self.obj_type, self.obj)
        
    def check_external_record_in_obj(self):
        '''Check if object (lane, library or sample) has external record for accession in self.accession'''
        ers = None
        if self.obj_type == 'lane':
            try:
                ers = ExternalRecord.objects.filter(lanes=self.obj, repository__name=self.repository_name)
            except ExternalRecord.DoesNotExist:
                LOGGER.info("No external records associated with lane %s.", self.obj)
        if self.obj_type == 'library':
            try:
                ers = ExternalRecord.objects.filter(libraries=self.obj, repository__name=self.repository_name)
            except ExternalRecord.DoesNotExist:
                LOGGER.info("No external records associated with library %s.", self.obj)
        if self.obj_type == 'sample':
            try:
                ers = ExternalRecord.objects.filter(samples=self.obj, repository__name=self.repository_name)
            except ExternalRecord.DoesNotExist:
                LOGGER.info("No external records associated with sample %s.", self.obj)

        found = False
        if len(ers) == 0:
            LOGGER.info("No ENA entries found!")
            return 0
        else:
            for er in ers:
                if er.accession == self.accession:
                    found = True
                    LOGGER.info("%s accession %s already associated with %s %s.", er.repository.name, er.accession, self.obj_type, self.obj)
                    return 1

        if not found:
            LOGGER.info("No entry for accession %s found!", self.accession)
            return 0

    def add_lane_obj(self, lane_id):
        try:
            self.obj = Lane.objects.get(id=lane_id)
        except Lane.DoesNotExist:
            LOGGER.error("No lane with id %d.", lane_id)
            sys.exit(1)
        self.obj_type = 'lane'

    def add_library_obj(self, code):
        try:
            self.obj = Library.objects.get(code=code)
        except Library.DoesNotExist:
            LOGGER.error("No library with code \"%s\".", code)
            sys.exit(1)
        self.obj_type = 'library'

    def add_sample_obj(self, sample_name):
        try:
            self.obj = Sample.objects.get(name=sample_name)
        except Sample.DoesNotExist:
            LOGGER.error("No sample with name \"%s\".", sample_name)
            sys.exit(1)
        self.obj_type = 'sample'

    def add_sample_obj_library(self, code):
        try:
            self.obj = Sample.objects.get(library__code=code)
        except Sample.DoesNotExist:
            LOGGER.error("No sample with name \"%s\".", sample_name)
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
            LOGGER.error("Failed to parse %s. Release_date is missing!", fname)
            failed = True
        if self.accession is None:
            LOGGER.error("Failed to parse %s. Acccession is missing!", fname)
            failed = True
        if self.receipt_type is None:
            LOGGER.error("Failed to parse %s. Receipt is of unknown type. (%s expected).", fname, " or ".join(receipt_types))
            failed = True
        if alias is None:
            LOGGER.error("Failed to parse %s. Failing to parse \"%s\" into anything meaningful.", fname, alias)
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
                LOGGER.error("No lane found for code=%s, facility=%s, lanenum=%d, flowcell=%s, flowlane=%s.", code, facility, lanenum, flowcell, flowlane)
                sys.exit(1)
            self.obj_type = 'lane'
        
        if self.receipt_type == 'EXPERIMENT':
            self.add_library_obj(code=alias)
            self.obj_type = 'library'
            
        if self.receipt_type == 'SAMPLE':
            self.add_sample_obj(sample_name=alias)
            self.obj_type = 'sample'
