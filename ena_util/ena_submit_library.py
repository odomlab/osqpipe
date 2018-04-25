#!/bin/env python
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
import datetime
from subprocess import Popen, PIPE
import xml.etree.ElementTree as ET

import django
from osqpipe.models import Lane, Lanefile, Library, SourceTreatment, ExternalRecord, ExternalRepository, Characteristic, Sample
from osqpipe.models import Machine

from logging import INFO, DEBUG
from osqutil.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)

django.setup()

class ExternalRecordManager(object):

  def __init__(self, entry_type=None, accession = None, release_date = None, is_public = False, repository_name='EBI ENA', code=None):

    self.receipt_type = None
    self.accession = None
    self.release_date = None
    self.obj_type = None
    self.obj = None
    self.repository_name = repository_name
    self.is_public = is_public
    self.code=code

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

#  Following has been deprecated as multiple samples may exist with same name but have different tissues associated. Use add_sample_obj_library(self, code) below instead!
#  def add_sample_obj(self, sample_name):
#    try:
#      self.obj = Sample.objects.get(name=sample_name)
#    except Sample.DoesNotExist:
#      LOGGER.error("No sample with name \"%s\".", sample_name)
#      sys.exit(1)
#    self.obj_type = 'sample'

  def add_sample_obj_library(self, code):
    try:
      self.obj = Sample.objects.get(library__code=code)
    except Sample.DoesNotExist:
      LOGGER.error("No sample with name \"%s\".", sample_name)
      sys.exit(1)
    self.obj_type = 'sample'

  def parse_xml_receipt(self, fname):
    '''Parse ENA submission receipt XML and adds accession from it to relevant lane/library/sample in repository'''

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
        # NB! External accessions are good but ENA own accessions are more universal for object modification and caceling. Hence 
        #     Following three lines are at the moment commented out!
        #for gchild in child:
        #  if gchild.tag == "EXT_ID":
        #    self.accession = gchild.get('accession')
      if child.tag == 'MESSAGES':
        messages = child
        errors = messages.findall('ERROR')
        if len(errors) > 0:
          accession_in_error = False
          for error in errors:
            LOGGER.error("%s\n", error.text)
            if "already exists as accession" in error.text:
              start = error.text.find('already exists as accession')
              self.accession =  error.text[(start+27):]
          if not accession_in_error:
            sys.exit("Invalid XML!\n\n")

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
    if self.receipt_type == 'SAMPLE':
      self.add_sample_obj_library(self.code)

class EnaAnnotation(object):
  '''A class coordinating extraction and formatting of meta-data necessary from ENA submission. This is done from lane perspective.'''
  
  def __init__(self, lane, external_repository, annotation_source='db', target_dir=""):

    self.target_dir = target_dir
    self.external_repository = external_repository
    
    # Define a list annotation values
    self.library_source = None
    self.library_selection = None
    self.library_strategy = None
    self.annotations = dict()
    self.taxon_id = None
    self.scientific_name = None
    self.common_name = None
    self.treatments_agent = []
    self.treatments_dose = []
    self.treatments_unit = []
    self.treatments_date = []
    self.grading = None
    self.code = None
    self.adapter1 = None
    self.adapter2 = None
    self.library_layout = None
    self.center_name = None
    self.run_center = None
    self.instrument_model = None
    self.facility = None
    self.facility_name = None
    self.flowlane = None
    self.lanenum = None
    self.rundate = None
    self.flowcell = None
    self.readlength = None
    self.filetype = None
    self.library_external_record = None
    self.sample_external_record = None
    self.files = []
    self.md5s = []

    # Define a secondary list of annotation values that are created through processing of the annotation values above
    self.study_alias = None
    self.sample_alias = None
    self.sample_title = None
    self.sample_description = None
    self.experiment_alias = None
    self.experiment_library_name = None
    self.experiment_library_construction_procotol = None
    self.experiment_title = None
    self.experiment_design_description = None
    self.run_alias = None
    # self.run_member_name = None

    self.library = None

    # Read annotations from source
    if annotation_source == 'db':
      self.get_annotations_from_db(lane=lane)
    else:
      sys.exit("Annotation source \"%s\" is not supported!\n\n")

    # Process primary annotation values to secondary required for ENA short read data submission
    self.process_annotations()

  def get_md5_file(self, fn):

    md5fn = fn + '.md5'
    if not os.path.isfile(md5fn):
      base = os.path.basename(fn)
      md5fn = os.path.join(self.target_dir, base + '.md5')
    
      if not os.path.isfile(md5fn):

        cmd = "md5sum %s > %s" % (fn, md5fn)

        LOGGER.info("Creating %s ...\n", md5fn)
        subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        (stdout, stderr) = subproc.communicate()
        retcode = subproc.wait()
        #if stdout is not None:
        #  sys.stdout.write(stdout)
        #if stderr is not None:
        #  sys.stderr.write(stderr)
        if retcode != 0:
          LOGGER.error("%s\nFailed to create %s.\n\n", stderr, md5fn)
          sys.exit(1)

    md5 = None
    LOGGER.info("Reading md5 sum from %s ...\n", md5fn)
    with open(md5fn,'rb') as fh:
      line = fh.readline().rstrip('\n')
      cols = line.split()
      md5 = cols[0]

    return md5

  def get_annotations_from_db(self, lane):
    '''Extract all annotations from db necessary to build an experiment (lane) submission'''

    self.library = Library.objects.get(id=lane.library_id)
    # Check if sample has been treated
    try:
      so = SourceTreatment.objects.get(source=self.library.sample.source)
    except SourceTreatment.DoesNotExist:
      so = None

    # Check that the library is either genome or exome
    if self.library.libtype.code not in ['genome', 'exome']:
      sys.exit("Library type '%s' not supported! (Lane_id=%d, library_code=%s)\n\n" % (self.library.libtype.code, lane.id, self.library.code))

    # Set values that are specific to whole genome sequencing and exome libraries
    self.library_source = 'GENOMIC'
    self.library_selection = 'RANDOM'
    if self.library.libtype.code == 'genome':
      self.library_strategy = 'WGS'
    if self.library.libtype.code == 'exome':
      self.library_strategy = 'WXS'
    self.filetype = 'fastq'
    self.sample_submitted = False

    # Check if generating sample XML is needed. This is done by checking if any other lanes of the associated library have already been
    # submitted to the repository. If so, the sample XML must have been submitted as well.
    if len(Lane.objects.filter(library=self.library.id, external_records__repository__name=self.external_repository)) > 0:
      self.sample_submitted = True
    # Sample annotations "ENA mutagenesis by carcinogen treatment checklist" 
    #    tissue_type
    #    sex
    #    date of birth
    #    date of death
    #    diagnosis
    #    strain
    #    tumor grading
    #    treatment agent
    #    treatment dose
    #    treatment date
    #    Further Details
    # Reference details related to a sample in form of an URI.

    self.annotations = dict()
    self.taxon_id = self.library.genome.species.accession
    self.scientific_name = self.library.genome.species.scientific_name
    self.common_name = self.library.genome.species.common_name
    if so is not None:
      self.treatments_agent.append(so.agent.name)
      self.treatments_dose.append(so.dose)
      if so.dose_unit is not None:
        self.treatments_unit.append(so.dose_unit.name)
      else:
        LOGGER.error("Dose unit for %s (sample %s) is missing. EXITING!\n" % (self.library.code, self.library.sample.name))
        sys.exit(1)
      self.treatments_date.append(so.date) # so.date.strftime("%Y-%m-%d_%H:%M:%S.%f")

    self.annotations['tissue_type'] = self.library.sample.tissue
    self.annotations['individual_name'] = self.library.sample.source.name
#   self.annotations['mother_name'] = self.library.sample.source.mother.name
#   self.annotations['father_name'] = self.library.sample.source.father.name
    self.annotations['strain'] = self.library.sample.source.strain.name
    self.annotations['date of birth'] = self.library.sample.source.date_of_birth
    self.annotations['date of death'] = self.library.sample.source.date_of_death
    self.annotations['sample name'] = self.library.sample.name
    if self.library.sample.source.sex is None:
      self.annotations['sex'] = 'NA'
    else:
      self.annotations['sex'] = self.library.sample.source.sex.name

    # Check if sample has diagnosis
    try:
      self.annotations['diagnosis'] = Characteristic.objects.get(category__iexact='Diagnosis',samples=self.library.sample).value
    except Characteristic.DoesNotExist:
      self.annotations['diagnosis'] = ''

    # Check if sample has TumourGrade
    try:
      self.annotations['tumour_grading'] = Characteristic.objects.get(category__iexact='TumourGrade',samples=self.library.sample).value
    except Characteristic.DoesNotExist:
      self.annotations['tumour_grading'] = ''

    # Check if sample has CausativeAgent
    try:
      self.annotations['causative_agent'] = Characteristic.objects.get(category__iexact='CausativeAgent',samples=self.library.sample).value
    except Characteristic.DoesNotExist:
      self.annotations['causative_agent'] = ''
    # Values from library and related
    self.code = self.library.code
    if self.library.adapter:
      self.adapter1 = self.library.adapter.sequence
    if self.library.adapter2:
      self.adapter2 = self.library.adapter2.sequence
    if lane.paired:
      self.library_layout = 'PAIRED'
    else:
      self.library_layout = 'SINGLE'
    # Values from lane and related
    # self.center_name = lane.facility.name
    self.center_name = 'CANCER RESEARCH UK CAMBRIDGE INSTITUTE' # NB! This is hard coded here as we are the data submitter!
    self.run_center = lane.facility.name # NB! This is not entirely correct! Correct run_center names need to be asked from ENA.
    self.instrument_model = lane.machine.platform
    self.facility = lane.facility.code
    self.facility_name = lane.facility.name
    self.flowlane = lane.flowlane
    self.lanenum = lane.lanenum
    self.rundate = lane.rundate
    self.flowcell = lane.flowcell
    self.readlength = lane.readlength
    for lanefile in Lanefile.objects.filter(lane=lane, filetype__name=self.filetype):
      fname = os.path.join(lanefile.archive.root_path, self.library.code, lanefile.filename + ".gz")
      if os.path.exists(fname):
        self.files.append(fname)
        # We can not use md5s from database as these are for uncompressed files. I.e md5 for gzipped file needs to be created from scratch.
        md5 = self.get_md5_file(fname)

        # self.md5s.append(lanefile.checksum)
        self.md5s.append(md5)
      else:
        sys.exit("File \"%s\" is not accessible or missing!\n\n" % fname)

  def process_annotations(self):
    '''Constructs values for missing ENA annotation categories.'''

    # Construct variables for ENA sample

    treatment_name = ''
    treatment_title = ''
    tumour_txt = ''
    if self.treatments_agent:
      treatment_name = '-'.join(self.treatments_agent)
      treatment_title = ' treated by ' + ' and treated by '.join(self.treatments_agent)
    else:
      treatment_name = 'noTreatment'

    diagnosis_txt = ''
    if self.annotations['diagnosis'] != '':
      treatment_txt = self.annotations['diagnosis']
      if self.annotations['diagnosis'] == 'HCC':
        treatment_txt = 'HCC nodule'
      if self.annotations['diagnosis'] == 'unclassified':
        treatment_txt = 'unclassified nodule'
    else:
      treatment_txt = 'sample'

    causative_txt = ' from individual'
    if self.annotations['causative_agent'] == 'spontaneous':
      causative_txt = ' of spontaneous origin from individual'
    if self.annotations['causative_agent'] == 'DEN':
      causative_txt = ' from DEN treated individual'

    # Construct sample alias as follows: donumber_tissue_treatment_individual_sciname.
    # Even though donumber is a value for library and in ENA/SRA vocabulary should belong to 'experiment', the donumber makes sample unique.
    # This is because in chipseq repository we do not distinguish between sample and library. Every library is assumed to have come from a unique library.
    # sample_alias = "%s_%s_%s_%s_%s" % (self.code, self.annotations['tissue'], treatment_name, self.annotations['sample name'], self.scientific_name)
    # self.sample_alias = sample_alias.replace(' ','_') 
    # self.sample_alias = self.annotations['individual_name']
    # self.sample_alias = self.annotations['sample name']
    self.sample_alias = "%s-%s" % (self.annotations['sample name'], self.annotations['tissue_type'])
    
    # Construct sample title
    self.sample_title = '%s %s %s%s %s' % (self.annotations['tissue_type'], treatment_txt, self.annotations['sample name'], causative_txt, self.annotations['individual_name'])
    # self.sample_title = "%s sample (%s) from individual %s%s" % (self.annotations['tissue_type'], self.annotations['sample name'], self.annotations['individual_name'], treatment_title)
    # Construct sample description. According to ENA/SRA, this should be a longer description of the sample, also showing how it differes from other samples.
    self.sample_description = self.sample_title + ". This sample is a unique specimen from this individual."
    
    # Construct variables for ENA experiment

    # construct run name from experiment name, facilitycode and flowlane
    # USE donumber instead of sample alias?
    # self.experiment_alias = self.sample_alias + "_" + self.library_layout + "_%s%d" % (self.facility, self.lanenum)
    # self.experiment_alias = "_".join([self.code, self.facility, "%d" % (self.lanenum)])
    self.experiment_alias = self.code
    
    # self.experiment_library_name = self.experiment_alias
    self.experiment_library_name = self.code
    if self.library_strategy == 'WGS':
      self.experiment_library_construction_procotol = 'DNA for the library preparation was extracted using QIAGEN AllPrep DNA/RNA mini kit using either single columns or 96 well plates. Library was constructed using Illumina TruSeq DNA PCR-Free Kit.'
    elif self.library_strategy == 'WXS':
      self.experiment_library_construction_procotol = 'Library was constructed using Illumina TruSeq DNA PCR-Free Kit.'
    else:
      sys.exit("Library strategy/type \"%s\" not supported!\n\n" % (self.library_strategy))

    # Build experiment title.
    self.experiment_title = "%s %sbp %s end sequencing of sample %s." % (self.instrument_model, str(self.readlength), self.library_layout.lower(), self.annotations['sample name'])
    
    # Build experiment description.
    adapters = ''
    if self.adapter1 is not None:
      adapters = 'index ' + self.adapter1 + '.'
    if self.adapter2 is not None:
      adapters = 'dual indexes ' + self.adapter1 + ' and ' + self.adapter2 + '.'
    self.experiment_design_description = self.experiment_title + ' Library was created with ' + adapters
    
    # Construct variables for ENA run
    self.run_alias = "_".join([self.code, self.facility, "%d" % (self.lanenum), "Flowcell=%s" % (self.flowcell), "Lane=%s" % str(self.flowlane)])
    # self.member_name = '%s_%d' % (self.flowcell, self.flowlane) # NB! member_name in DATA_BLOCK is not widely used. Neither does it show the multiplexing. ENA people say that it shows SOMETHING ELSE.

    # Find external records for library (experiment) and sample 
    try:
      ers = ExternalRecord.objects.get(libraries__code=self.experiment_alias, repository__name=self.external_repository)
      self.library_external_record = ers.accession
    except ExternalRecord.DoesNotExist:
      LOGGER.info("No external records for library \'%s\'.", self.experiment_alias)
    try:
      ers = ExternalRecord.objects.get(samples__id=self.library.sample.id, repository__name=self.external_repository)
      self.sample_external_record = ers.accession
    except ExternalRecord.DoesNotExist:
      LOGGER.info("No external records for sample \'%s\'.", self.sample_alias)      

class EnaXmlObject(object):
  '''A class / container to keep info about generated XML file, its type, alias, validation, submission etc. state'''
  
  def __init__(self, otype, alias, path=None, template_path=None, validated=False, submitted=False, modify=False):

    if otype not in ['experiment','sample','run']:
      LOGGER.error("Unknown ENA xml submission type \'%s\'. Exiting!", otype)
    
    self.otype = otype # object type. One of following values: experiment, sample, run
    self.alias = alias # experiment, sample or run alias
    self.validated = validated # logical. Shows if XML has been validated
    self.submitted = submitted # logical. Shows if XML has been submitted
    self.accession = None # ENA accession
    self.modify = modify
    fixed_alias = self.fix_alias_for_fname(alias) # alias that has some characters replaced so that it can be used in file names
    
    # check if XML template exists
    self.template = os.path.join(template_path, otype + ".xml")
    if not os.path.exists(self.template):
      LOGGER.error("XML template file %s missig or not accessible!\n\n", self.template)
      sys.exit(1)

    # check if submission XML template exists.
    submission_template_fn = 'submission.xml'
    if self.modify:
      submission_template_fn = 'submission_update.xml'
    self.submission_template = os.path.join(template_path, submission_template_fn)
    if not os.path.exists(self.submission_template):
      LOGGER.error("XML template file %s missig or not accessible!\n\n", self.submission_template)
      sys.exit(1)    

    # build XML and submission XML file names
    self.filename = os.path.join(path, fixed_alias + '-' + os.path.basename(self.template))
    self.submission_filename = os.path.join(path, fixed_alias + "-submission.xml")

    # create variables to keep data files and md5sums when otype='run'
    self.datafiles = []
    self.datamd5s = []
    
  def fix_alias_for_fname(self, alias):
    '''Fixes unwanted characters such as spaces and forward slashes in alias so that alias could be used in constructing filename'''
    alias = alias.replace(" ","_")
    alias = alias.replace("/","_")
    return alias
    
class EnaXmlCreator(object):
  '''A class coordinating construction of XML documents and submission XML files for ENA.'''
  
  def __init__(self, lane, study_alias, template_dir="", target_dir="", release_date_str=None, repository_name = 'EBI ENA', modify=False):
    '''Takes study alias and optionally, path to xml temlates (template_dir) and targets (dir where xml files should be written).''' 

    self.study_alias = study_alias
    self.template_dir = template_dir
    self.target_dir = target_dir
    self.release_date_str = None
    self.external_repository = repository_name
    self.modify = modify # if set, the submission XML will be for submission of modification/replacement XML. I.e. instead of ADD, MODIFY will be used.

    # Set release date in case not provided, assign it 2 years from now
    if release_date_str is None:
      today = datetime.date.today()
      self.release_date_str = today.replace(year=today.year + 2)
    else:
      self.release_date_str = datetime.strptime(ARGS.release_date , '%Y-%m-%d')
    
    # read annotations from repository
    self.a = EnaAnnotation(lane, self.external_repository, target_dir=self.target_dir)
    self.a.process_annotations()
    self.a.study_alias = study_alias

  def write_experiment_xml(self):
    '''Writes experiment xml file for library. Returns ENA xml info object.'''

    # Build experiment xml info object
    exo = EnaXmlObject('experiment', self.a.experiment_alias, path=self.target_dir, template_path=self.template_dir, modify=self.modify)

    if self.a.library_external_record is not None:
      exo.submitted = True
      exo.accession = self.a.library_external_record
      LOGGER.info("Experiment already submitted (accession=%s)", exo.accession)
      if not self.modify:
        return exo
        
    # Check if XML may have already been created
    if os.path.isfile(exo.filename):
      if self.modify:
        LOGGER.info("XML for %s \'%s\' already exists: %s. Owerwriting XML!", exo.otype, exo.alias, exo.filename)
        # sys.exit(1)
      else:
        LOGGER.info("XML for %s \'%s\' already exists: %s", exo.otype, exo.alias, exo.filename)
      return exo
    
    # parse template
    xml = ET.ElementTree()
    xml.parse(exo.template)

    # Modify values in template
    experiment = xml.find('EXPERIMENT')
    experiment.set('alias',exo.alias)
    experiment.set('center_name',self.a.center_name)
    title = experiment.find('TITLE')
    title.text = self.a.experiment_title
    # print "STUDY_REF=%s" % self.a.study_alias
    experiment.find('STUDY_REF').set('refname', self.a.study_alias)
    design = experiment.find('DESIGN')
    design.find('DESIGN_DESCRIPTION').text = self.a.experiment_design_description
    design.find('SAMPLE_DESCRIPTOR').set('refname', self.a.sample_alias)
    ldesc = design.find('LIBRARY_DESCRIPTOR')
    ldesc.find('LIBRARY_NAME').text = self.a.experiment_library_name
    ldesc.find('LIBRARY_STRATEGY').text = self.a.library_strategy
    ldesc.find('LIBRARY_SOURCE').text = self.a.library_source
    ldesc.find('LIBRARY_SELECTION').text = self.a.library_selection
    library_layout = ldesc.find('LIBRARY_LAYOUT')
    # remove all existing elements under the library layout
    for child in library_layout:
      library_layout.remove(child)
    # add one and only element:
    library_layout.append( ET.Element(self.a.library_layout) )
    ldesc.find('LIBRARY_CONSTRUCTION_PROTOCOL').text = self.a.experiment_library_construction_procotol
    platform = experiment.find('PLATFORM')
    illumina = platform.find('ILLUMINA')
    illumina.find('INSTRUMENT_MODEL').text = self.a.instrument_model
    
    # write xml into file
    xml.write(exo.filename)

    # write submission xml
    self.write_submission_xml(exo)

    return exo
  
  def write_run_xml(self):
    '''Reads ENA run XML from xml_template and writes run XML for the lane. Returns ENA xml info object.'''

    # Build run xml info object
    exo = EnaXmlObject('run', self.a.run_alias, path=self.target_dir, template_path=self.template_dir, modify=self.modify)
    
    # Check if XML may have already been created
    if os.path.isfile(exo.filename):
      LOGGER.info("XML for %s \'%s\' already exists: %s", exo.otype, exo.alias, exo.filename)
      # Add filenames to list of datafiles
      for i in range(0,len(self.a.files)):
        exo.datafiles.append(self.a.files[i])
        exo.datamd5s.append(self.a.md5s[i])

      if not self.modify:
        return exo

    # Read XML template
    xml = ET.ElementTree()
    xml.parse(exo.template)

    # edit xml template
    run = xml.find('RUN')
    # set run name/alias
    run.set('alias', self.a.run_alias)
    run.set('center_name', self.a.center_name)
    run.set('run_center', self.a.center_name)
    run.set('run_date', str(self.a.rundate) + 'T10:00:00')
    run.find('EXPERIMENT_REF').set('refname', self.a.experiment_alias)
    data_block = run.find('DATA_BLOCK')    
    # data_block.set('member_name', self.a.member_name) # According to ENA staff, this is not used for multiplexing.
    # member_name="TODO: FOR DEMULTIPLEXED DATA ONLY (see note below)"
    data_block.set('member_name','')

    files = data_block.find('FILES')
    f = files.find('FILE')
    f.set('filename', os.path.basename(self.a.files[0]))
    f.set('filetype', self.a.filetype)
    f.set('checksum', self.a.md5s[0])
    exo.datafiles.append(self.a.files[0])
    exo.datamd5s.append(self.a.md5s[0])
    for i in range(1,len(self.a.files)):
      f = ET.Element('FILE', attrib={ 'filename':os.path.basename(self.a.files[i]), 'checksum_method':"MD5", 'filetype':self.a.filetype, 'checksum':self.a.md5s[i]})
      files.append(f)
      exo.datafiles.append(self.a.files[i])
      exo.datamd5s.append(self.a.md5s[i])

    # write XML
    xml.write(exo.filename)

    # write submission xml
    self.write_submission_xml(exo)
    
    return exo
    
  def write_sample_xml(self):
    '''Reads ENA sample XML from template and writes sample XML. Returns ENA xml info object.'''

    # Build run xml info object
    exo = EnaXmlObject('sample', self.a.sample_alias, path=self.target_dir, template_path=self.template_dir, modify=self.modify)

    if self.a.sample_external_record is not None:
      exo.submitted = True
      exo.accession = self.a.sample_external_record

      LOGGER.info("Sample \'%s\' already submitted under accession \'%s\'", self.a.sample_alias, self.a.sample_external_record)
      if not self.modify:
        return exo

    # Check if XML may have already been created
    if os.path.isfile(exo.filename):
      if self.modify:
        LOGGER.info("XML for %s \'%s\' already exists: %s. Remove the file before XML update could be written!", exo.otype, exo.alias, exo.filename)
        sys.exit(1)
      else:
        LOGGER.info("XML for %s \'%s\' already exists: %s", exo.otype, exo.alias, exo.filename)
      return exo

    # Read XML template
    xml = ET.ElementTree()
    xml.parse(exo.template)
    
    # edit xml template
    sample = xml.find('SAMPLE')
    # set sample name    
    sample.set('alias', self.a.sample_alias)
    sample.set('center_name', self.a.center_name)
    # set sample title
    sample.find('TITLE').text = self.a.sample_title
    # set sample taxon, scientific name and common name
    sample_name = sample.find('SAMPLE_NAME')
    sample_name.find('TAXON_ID').text = str(self.a.taxon_id)
    sample_name.find('SCIENTIFIC_NAME').text = self.a.scientific_name
    sample_name.find('COMMON_NAME').text = self.a.common_name
    # set sample description
    sample.find('DESCRIPTION').text = self.a.sample_description
    # add sample annotations
    sample_attributes = sample.find('SAMPLE_ATTRIBUTES')
    # Identify sample_attribute. The code below assumes only one SAMPLE_ATTRIBUTE being defined in the template
    # and self.annotations containing at least one category-value pair.    
    sample_attribute_list = sample_attributes.findall('SAMPLE_ATTRIBUTE')    
    # Remove all existing sample attributes
    for sa in sample_attribute_list:
      sample_attributes.remove(sa)
    for a in self.a.annotations:
      if a == 'causative_agent':
        continue
      sample_attribute = ET.Element('SAMPLE_ATTRIBUTE')
      ET.SubElement(sample_attribute, 'TAG').text = a
      ET.SubElement(sample_attribute, 'VALUE').text = str(self.a.annotations[a])
      sample_attributes.append(sample_attribute)
    if len(self.a.treatments_agent) > 0:
      for agent,dose,unit,date in zip(self.a.treatments_agent, self.a.treatments_dose, self.a.treatments_unit, self.a.treatments_date):
        sample_attribute = ET.Element('SAMPLE_ATTRIBUTE')
        ET.SubElement(sample_attribute, 'TAG').text = 'treatment agent'
        ET.SubElement(sample_attribute, 'VALUE').text = agent
        sample_attributes.append(sample_attribute)
        sample_attribute = ET.Element('SAMPLE_ATTRIBUTE')
        ET.SubElement(sample_attribute, 'TAG').text = 'treatment dose'
        ET.SubElement(sample_attribute, 'VALUE').text = str(dose)
        ET.SubElement(sample_attribute, 'UNITS').text = unit
        sample_attributes.append(sample_attribute)        
        sample_attribute = ET.Element('SAMPLE_ATTRIBUTE')
        ET.SubElement(sample_attribute, 'TAG').text = 'treatment date'
        ET.SubElement(sample_attribute, 'VALUE').text = str(date)
        sample_attributes.append(sample_attribute)
    else:
      sample_attribute = ET.Element('SAMPLE_ATTRIBUTE')
      ET.SubElement(sample_attribute, 'TAG').text = 'treatment agent'
      ET.SubElement(sample_attribute, 'VALUE').text = 'not treated'
      sample_attributes.append(sample_attribute)

    # write xml into file
    xml.write(exo.filename)

    # write submission xml
    self.write_submission_xml(exo)

    return exo

  def write_submission_xml(self, exo):
    '''Takes ENA xml info object and creates submission XML.'''

    # Check if submission XML may have already been created
    if os.path.isfile(exo.submission_filename):
      LOGGER.info("Submission XML for %s \'%s\' already exists: %s", exo.otype, exo.alias, exo.submission_filename)
      return exo
    
    # Read submission XML template
    xml = ET.ElementTree()
    xml.parse(exo.submission_template)

    action_type = 'ADD'
    if self.modify:
      action_type = 'MODIFY'
    
    submission = xml.find('SUBMISSION')
    submission.set('alias', exo.alias + "-" + exo.otype + "_submission")        
    submission.set('center_name',self.a.center_name)
    actions = submission.find('ACTIONS')
    i = 1
    for action in actions.findall('ACTION'):
      for child in action:
        if child.tag == action_type:
          if child.attrib['schema'] == exo.otype:
            child.set('source', exo.filename)
          else:
            #if child.attrib['schema'] not in ['run','experiment','sample']:
            # sys.stderr.write('Unsupported ACTION element with schema=\'%s\' in template! Removing the corresponding ACTION element.\n' % child.attrib['schema'])
            actions.remove(action)
        if child.tag == 'HOLD':        
          child.set('HoldUntilDate', self.release_date_str.strftime("%Y-%m-%d"))        

    # write xml into file
    xml.write(exo.submission_filename)

class enaXMLuploader(object):
  '''A class coordinating submission of XML documents to ENA'''
  def __init__(self, credentials_file):

    self.credentials = None
    self.read_credentials(credentials_file)

  def read_credentials(self, credentials_file):
    '''Parse ENA credentials file.'''

    self.credentials = {'username':None, 'password':None}
    with open(credentials_file,'rb') as fh:
      for line in fh:
        line = line.rstrip('\n')
        cols = line.split('\t')
        if cols[0] in self.credentials:
          self.credentials[cols[0]] = cols[1]
      fh.close()
    for key in self.credentials:
      if self.credentials[key] is None:
        LOGGER.error("Ill formated credentials file. No value for \'%s\'!\n\n", key)
        sys.exit(1)

  def upload_xml(self, s_xml, f_xml=None, s_type=None, validate=True):
    '''Submits XML to ENA, parses result for failure or of accession. With validation_only being set, submission is directed instead to ENA test server fro XML syntax validation.'''
    # s_xml - submission XML for f_xml file.
    # f_xml - XML file to be submitted.
    # s_type - submission type e.g. ['RUN','ANALYSIS','SAMPLE','EXPERIMENT']
    # validate - if set the submission is directed to ENA test server for XML syntax validation only.
    
    r_fname = s_xml + ".receipt"

    LOGGER.info('Preparing upload of following XMLs to ENA using curl.')
    cmd = "curl -k"
    cmd += " -F \"SUBMISSION=@" + s_xml + "\" "
    if s_type is not None and f_xml is not None:
      cmd += "-F \"" + s_type.upper() + "=@" + f_xml + "\""
    if validate:
      cmd += " \"https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%20" + self.credentials['username'] + "%20" + self.credentials['password'] + "\""
      r_fname = s_xml + ".test-receipt"
    else:
      cmd += " \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%20" + self.credentials['username'] + "%20" + self.credentials['password'] + "\""
    LOGGER.info(cmd)
    if os.path.isfile(r_fname):
      LOGGER.info("Skipping, already submitted.")
      if os.path.isfile(r_fname):
        LOGGER.info("Found submission receipt %s.", r_fname)
        return r_fname
      else:
        LOGGER.error("%s expected but not found!\n\n", r_fname)
        return ""
    subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    (stdout, stderr) = subproc.communicate()
    retcode = subproc.wait()
    if stderr is not None:
      LOGGER.info("STDERR:\n%s\n", stderr)
    if retcode != 0:
      LOGGER.error("%sXML transfer failed!\n\n", stderr)
      sys.exit(1)
    if stdout is not None:
      # Write receipt on disk
      with open(r_fname, "wb") as fout:
        fout.write(stdout)
        fout.close()
      # Write receipt to stdout as well
      LOGGER.info("STDOUT:\n%s\n", stdout)

      # return receipt file name
      return r_fname
    else:
      sys.exit("No receipt received for command \"%s\"" % cmd)

  def upload_cancel_xml(self, s_xml, validate=True):
    '''Submits object cancellation XML to ENA'''
    # TODO: currently receipt file name is not catched and analysed for the result!
    self.upload_xml(self, s_xml, validate=True)

  def submit_files(self, file_dict, s_type, test=True):
    '''Submits all files in file dict. File dict contains file names as keys and sublission XML file names for the files as values'''
    for f in file_dict:
      receipt_file = self.upload_xml(file_dict[f], f, s_type, validate=test)
      LOGGER.info("Receipt=%s" % receipt_file)

  def submit_files_in_directory(self, path, xml_type, validate_only):
    '''Reads files in a directory, identifies ENA submission files by file suffixes and submit these to ENA in a specific order.'''
    
    os.chdir(path)
    files = os.listdir(path)

    samples = dict()
    experiments = dict()
    runs = dict()
    analysis = dict()

    # Read files in a directory and sort them by file type to 4 dictionaries
    for f in files:
      if f.endswith('-sample.xml'):
        s = f.replace('-sample.xml','-submission.xml')
        samples[f] = s
        if not os.path.isfile(os.path.join(path,s)):
          LOGGER.error("%s missing!\n", s)
          sys.exit(1)
      if f.endswith('-experiment.xml'):
        s = f.replace('-experiment.xml','-submission.xml')
        experiments[f] = s
        if not os.path.isfile(os.path.join(path,s)):
          LOGGER.error("%s missing!\n", s)
          sys.exit(1)
      if f.endswith('-run.xml'):
        s = f.replace('-run.xml','-submission.xml')
        runs[f] = s
        if not os.path.isfile(os.path.join(path,s)):
          LOGGER.error("%s missing!\n", s)
          sys.exit(1)
      if f.endswith('-analysis.xml'):
        s = f.replace('-analysis.xml','-submission.xml')
        analysis[f] = s
        if not os.path.isfile(os.path.join(path,s)):
          LOGGER.error("%s missing!\n", s)
          sys.exit(1)

      if xml_type == "SAMPLE":
        self.submit_files(samples, "SAMPLE", test=validate_only)
      if xml_type == "EXPERIMENT":
        self.submit_files(experiments, "EXPERIMENT", test=validate_only)
      if xml_type == "RUN":
        self.submit_files(runs, "RUN", test=validate_only)
      if xml_type == "ANALYSIS":
        self.submit_files(runs, "ANALYSIS", test=validate_only)

class enaFileUploader(object):
  '''A class for uploading data files to ENA Aspera server'''
  # NB! In future it would be good to abstract this class to general Aspera communication.
  
  def __init__(self, credentials_file, files, submission_dir):
    '''Computes md5 sums and uploads files with sums to ENA webin using credentials from credentials_file.'''

    credentials = self.read_credentials(credentials_file)
    os.environ['ASPERA_SCP_PASS'] = credentials['password']

    for fn in files:
      
      #md5fn = ""
      #base = os.path.basename(fn)
      #if submission_dir is None:
      #  temp_dir = tempfile.gettempdir()
      #  md5fn = os.path.join(temp_dir, base + '.md5')
      #else:
      #  md5fn = os.path.join(submission_dir, base + '.md5')
      #  md5f = create_md5_file(fn, md5fn)

      #self.ena_upload_file(credentials, md5f)
      # os.unlink(md5f)                                                                                                                                                           
      self.ena_upload_file(credentials, fn)
      
  def read_credentials(self, credentials_file):
    '''Parse ENA credentials file.'''

    credentials = {'username':None, 'password':None}
    with open(credentials_file,'rb') as fh:
      for line in fh:
        line = line.rstrip('\n')
        cols = line.split('\t')
        if cols[0] in credentials:
          credentials[cols[0]] = cols[1]
      fh.close()
    for key in credentials:
      if credentials[key] is None:
        LOGGER.error("Ill formated credentials file. No value for \'%s\'!\n\n", key)
        sys.exit(1)

    return credentials

  def create_md5_file(self, fn, md5fn):
    '''Creates md5 for a file'''

    if os.path.isfile(fn):
      LOGGER.error("%s already exists.\n", md5fn)
      return fn

    cmd = "md5sum %s > %s" % (fn, md5fn)

    LOGGER.info("Creating %s ...\n", md5fn)
    subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    (stdout, stderr) = subproc.communicate()
    retcode = subproc.wait()

    if retcode != 0:
      LOGGER.error("%s\nFailed to create %s.\n\n", stderr, md5fn)
      sys.exit(1)

    return md5fn

  def ena_upload_file(self, credentials, fn):
    '''Uses ascp to upload a single data file to ENA webin using credentials provided.'''
    
    retry = 3
    
    while(retry > 0):
      # Note that we are uploading gzipped file.
      cmd = "ascp -QT -l300M -L- " + fn + " " + credentials['username'] + "@webin.ebi.ac.uk:."
      LOGGER.info("Uploading %s ...\n" % fn)
      LOGGER.info(cmd)
      subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
      (stdout, stderr) = subproc.communicate()
      retcode = subproc.wait()
      
      if retcode != 0:
        if "Data transfer stalled, timed out" in stderr:
          LOGGER.error("File transfer stalled! Re-starting file transfer.\n\n")
          retry = retry -1
        else:          
          LOGGER.error("%s\nFile transfer failed!\n\n", stderr)
          sys.exit(1)
      else:
        return
      
class EnaSubmitter(object):
  '''A class coordinating data submission to ENA by: triggering XML file generation, uploading data, uploading XMLs, parsing submission accession and entering these to repository'''

  def __init__(self, study_alias, study_xml, xml_template_dir, credentials_file, release_date_str=None, validate_xml=True, upload_files=False, submit_xml=False, target_dir="", force=False, modify=False):
    '''Takes some variables to set up the environment for the submission.'''

    self.external_repository = 'EBI ENA'
    
    self.study_alias = study_alias
    self.study_xml = study_xml
    self.xml_template_dir = xml_template_dir
    self.credentials_file = credentials_file
    self.release_date_str = release_date_str
    self.validate_xml = validate_xml
    self.upload_files = upload_files
    self.submit_xml = submit_xml
    self.target_dir = target_dir
    self.force = force
    self.modify = modify
    
  def submit_lane(self, lane, create_sample_xml=True, create_experiment_xml=True, create_run_xml=True):
    '''Submits data and meta-data for a particular lane_id if not submitted before.'''

    # install XML creator
    ex = EnaXmlCreator(lane, self.study_alias, self.xml_template_dir, self.target_dir, repository_name = self.external_repository, modify=self.modify)

    # create XMLs
    if create_sample_xml:
      sample_xml = ex.write_sample_xml()
    if create_experiment_xml:
      experiment_xml = ex.write_experiment_xml()
    if create_run_xml:
      run_xml = ex.write_run_xml()

    # Upload data files
    if self.upload_files:
      for datafile in run_xml.datafiles:
        efu = enaFileUploader(self.credentials_file, [datafile], os.path.split(datafile)[0])

    if not self.validate_xml and not self.submit_xml:
      return
    # Install XML uploader
    ena_loader = enaXMLuploader(self.credentials_file)

    # If XMLs have not been submitted then validate and if necessary submit as well.
    if create_sample_xml:
      if not sample_xml.submitted or self.modify:
        # validate XML
        if self.validate_xml:
          receipt_file = ena_loader.upload_xml(sample_xml.submission_filename, sample_xml.filename, 'SAMPLE', validate=True)
        # submit XML
        if self.submit_xml:
          receipt_file = ena_loader.upload_xml(sample_xml.submission_filename, sample_xml.filename, 'SAMPLE', validate=False)
          erm = ExternalRecordManager(sample_xml.otype, release_date = self.release_date_str, is_public = False, repository_name = self.external_repository, code=ex.a.code)
          # Extract assigned accession from receipt and insert it to repository
          erm.parse_xml_receipt(receipt_file)
          erm.add_external_record_to_obj()

    if create_experiment_xml:
      if not experiment_xml.submitted or self.modify:
        # validate XML
        if self.validate_xml:
          receipt_file = ena_loader.upload_xml(experiment_xml.submission_filename, experiment_xml.filename, 'EXPERIMENT', validate=True)
        # submit XML
        if self.submit_xml:
          receipt_file = ena_loader.upload_xml(experiment_xml.submission_filename, experiment_xml.filename, 'EXPERIMENT', validate=False)
          erm = ExternalRecordManager(experiment_xml.otype, release_date = self.release_date_str, is_public = False, repository_name = self.external_repository)
          # Extract assigned accession from receipt and insert it to repository
          erm.parse_xml_receipt(receipt_file)
          erm.add_external_record_to_obj()
    if create_run_xml:
      if not run_xml.submitted or self.modify:
        # validate XML
        if self.validate_xml:
          receipt_file = ena_loader.upload_xml(run_xml.submission_filename, run_xml.filename, 'RUN', validate=True)
        # submit XML
        if self.submit_xml:
          receipt_file = ena_loader.upload_xml(run_xml.submission_filename, run_xml.filename, 'RUN', validate=False)
          erm = ExternalRecordManager(run_xml.otype, release_date = self.release_date_str, is_public = False, repository_name = self.external_repository)
          # Extract assigned accession from receipt and insert it to repository
          erm.parse_xml_receipt(receipt_file)
          erm.add_external_record_to_obj()

  def submit_project(self, project):
    '''Submits data and meta-data for all lanes not yet been submitted for a project.'''

    # Find lanes of a project with no association to external repository 'EBI ENA'  
    try:
      lanes = Lane.objects.filter(library__projects__name=project).exclude(external_records__repository__name=self.external_repository)
    except Lane.DoesNotExist:
      LOGGER.error("No lanes for library code=%d!", libcode)
      sys.exit(1)

    # Submit data for one lane at the time
    for lane in lanes:
      self.submit_lane(lane)
        
  def submit_lane_with_id(self, lane_id):
    '''Submits data and meta-data for lane given its id.'''

    try:
      lane = Lane.objects.get(id=lane_id)
    except Lane.DoesNotExist:
      LOGGER.error("No lane with id=%d! Exiting!", lane_id)
      sys.exit(1)
    self.submit_lane(lane)
        
  def submit_library(self, libcode):
    '''Submits data and meta-data for library'''

    try:
      lanes = Lane.objects.filter(library__code=libcode).exclude(external_records__repository__name=self.external_repository)
    except Lane.DoesNotExist:
      LOGGER.error("No lanes for library code=%d!", libcode)
      sys.exit(1)    

    for lane in lanes:
      self.submit_lane(lane)

  def create_sample_xml(self, libcode):
    '''Creates and submits sample XML only.'''

    try:
      if self.force or self.modify:
        print "Forcing database query by ignoring associations to existing external records."
        lanes = Lane.objects.filter(library__code=libcode)
      else:
        lanes = Lane.objects.filter(library__code=libcode).exclude(external_records__repository__name=self.external_repository)
    except Lane.DoesNotExist:
      LOGGER.error("No lanes for library code=%d!", libcode)
      sys.exit(1)
      
    if len(lanes) > 0:
      self.submit_lane(lanes[0], create_experiment_xml=False, create_run_xml=False)
    else:
      LOGGER.error("Number of lanes retrieved was 0")
      sys.exit(1)

  def create_experiment_xml(self, libcode):
    '''Creates and submits experiment XML only.'''

    try:
      if self.force or self.modify:
        print "Forcing database query by ignoring associations to existing external records."
        lanes = Lane.objects.filter(library__code=libcode)
      else:
        lanes = Lane.objects.filter(library__code=libcode).exclude(external_records__repository__name=self.external_repository)
    except Lane.DoesNotExist:
      LOGGER.error("No lanes for library code=%d!", libcode)
      sys.exit(1)
    
    if len(lanes) > 0:
      self.submit_lane(lanes[0], create_sample_xml=False, create_run_xml=False)
    else:
      LOGGER.error("Number of lanes retrieved was 0")
      sys.exit(1)

  def submit_lanes_from_file(self, fn):
    '''Submits data and meta-data for a set of lane_ids found in the first column of the file fn'''
    with open(fn, 'rb') as fh:
      for line in fh:
        line = line.rstrip('\n')
        cols = line.split('\t')
        self.submit_lane_with_id(cols[0])

if __name__ == '__main__':

  from argparse import ArgumentParser

  # Required arguments
  
  PARSER = ArgumentParser(description='Data submission to pre-exisgting ENA study.')

  # Define what is being submitted:
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

  # Submit sample only
  PARSER.add_argument('--create_sample_xml', dest='create_sample_xml', action='store_true', default=False,
                      help='Create sample XML only.')
  PARSER.add_argument('--create_experiment_xml', dest='create_experiment_xml', action='store_true', default=False,
                      help='Create experiment XML only.')
  PARSER.add_argument('--force', dest='force', action='store_true', default=False,
                      help='Forces file generation even if accession in repository already exists.')
  PARSER.add_argument('--modify', dest='modify', action='store_true', default=False,
                      help='XML is created for an update of object already in ENA.')

  # Define ENA study where the data should be added. This can be done either by:
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
  # Define the action required from the tool:
  # - upload all files (and create md5 sums if not present in the file path)
  # - create XMLs
  # - validate XMLs
  # - submit XMLs
  #PARSER_ACTION_GROUP = PARSER.add_mutually_exclusive_group(required=True)
  #PARSER_ACTION_GROUP.add_argument('--upload_files', dest='upload_files', action='store_true',
  #                    help='Upload files for submission')  
  #PARSER_ACTION_GROUP.add_argument('--create_xml', dest='create_xml', action='store_true',
  #                    help='Create and validate XMls for submission.')
  #PARSER_ACTION_GROUP.add_argument('--submit_xml', dest='submit_xml', action='store_true',
  #                    help='Submit XMls to ENA.')
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
  ena = EnaSubmitter(study_alias=ARGS.study_alias, study_xml=ARGS.study_xml, xml_template_dir=ARGS.xml_template_dir, credentials_file=ARGS.ena_credentials_file, release_date_str=ARGS.release_date, validate_xml=ARGS.validate_xml, upload_files=ARGS.upload_files, submit_xml=ARGS.submit_xml, target_dir=ARGS.target_dir, force=ARGS.force, modify=ARGS.modify)
  
  if ARGS.project:
    ena.submit_project(project=ARGS.project, )

  if ARGS.lane_id:
    # e.g. lane_id=75696
    # e.g. lane_id=77516
    ena.submit_lane_with_id(ARGS.lane_id)
  
  if ARGS.lane_id_file:
    ena.submit_lanes_from_file(fn=ARGS.lane_id_file)

  if ARGS.libcode:
    if ARGS.create_sample_xml:
      ena.create_sample_xml(libcode=ARGS.libcode)
    elif ARGS.create_experiment_xml:
      ena.create_experiment_xml(libcode=ARGS.libcode)
    else:
      ena.submit_library(libcode=ARGS.libcode)    

## Upload associated files
## Generate and validate XML files (in target directory)
## Upload XML files (in target directory)
