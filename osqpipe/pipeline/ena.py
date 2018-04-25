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
import datetime
from subprocess import Popen, PIPE
import xml.etree.ElementTree as ET

import django
from osqpipe.models import Lane, Lanefile, Library, SourceTreatment, \
  ExternalRecord, ExternalRepository, Characteristic, Sample, Machine
from osqpipe.pipeline.external_record import ExternalRecordManager

from logging import INFO, DEBUG
from osqutil.setup_logs import configure_logging
LOGGER = configure_logging(level=INFO)

django.setup()

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

    # Read annotations from source
    if annotation_source == 'db':
      self.get_annotations_from_db(lane=lane)
    else:
      sys.exit('Annotation source \'%s\' is not supported!\n\n')

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
    with open(md5fn,'rb') as fh:
      line = fh.readline().rstrip('\n')
      cols = line.split()
      md5 = cols[0]

    return md5

  def get_annotations_from_db(self, lane):
    '''Extract all annotations from db necessary to build an experiment (lane) submission'''

    library = Library.objects.get(id=lane.library_id)
    # Check if sample has been treated
    try:
      so = SourceTreatment.objects.get(source=library.sample.source)
    except SourceTreatment.DoesNotExist:
      so = None

    # Check that the library is either genome or exome
    if library.libtype.code not in ['genome', 'exome']:
      sys.exit("Library type '%s' not supported! (Lane_id=%d, library_code=%s)\n\n" % (library.libtype.code, lane.id, library.code))

    # Set values that are specific to whole genome sequencing and exome libraries
    self.library_source = 'GENOMIC'
    self.library_selection = 'RANDOM'
    if library.libtype.code == 'genome':
      self.library_strategy = 'WGS'
    if library.libtype.code == 'exome':
      self.library_strategy = 'WXS'
    self.filetype = 'fastq'
    self.sample_submitted = False

    # Check if generating sample XML is needed. This is done by checking if any other lanes of the associated library have already been
    # submitted to the repository. If so, the sample XML must have been submitted as well.
    if len(Lane.objects.filter(library=library.id, external_records__repository__name=self.external_repository)) > 0:
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
    self.taxon_id = library.genome.species.accession
    self.scientific_name = library.genome.species.scientific_name
    self.common_name = library.genome.species.common_name
    if so:
      self.treatments_agent.append(so.agent.name)
      self.treatments_dose.append(so.dose)
      self.treatments_unit.append(so.dose_unit.name)
      self.treatments_date.append(so.date) # so.date.strftime("%Y-%m-%d_%H:%M:%S.%f")

    self.annotations['tissue_type'] = library.sample.tissue
    self.annotations['individual_name'] = library.sample.source.name
#   self.annotations['mother_name'] = library.sample.source.mother.name
#   self.annotations['father_name'] = library.sample.source.father.name
    self.annotations['strain'] = library.sample.source.strain.name
    self.annotations['date of birth'] = library.sample.source.date_of_birth
    self.annotations['date of death'] = library.sample.source.date_of_death
    self.annotations['sample name'] = library.sample.name
    self.annotations['sex'] = library.sample.source.sex.name

    # Check if sample has diagnosis
    try:
      self.annotations['diagnosis'] = Characteristic.objects.get(category__iexact='Diagnosis',samples=library.sample).value
    except Characteristic.DoesNotExist:
      self.annotations['diagnosis'] = ''

    # Check if sample has TumourGrade
    try:
      self.annotations['tumour_grading'] = Characteristic.objects.get(category__iexact='TumourGrade',samples=library.sample).value
    except Characteristic.DoesNotExist:
      self.annotations['tumour_grading'] = ''

    # Check if sample has CausativeAgent
    try:
      self.annotations['causative_agent'] = Characteristic.objects.get(category__iexact='CausativeAgent',samples=library.sample).value
    except Characteristic.DoesNotExist:
      self.annotations['causative_agent'] = ''
    # Values from library and related
    self.code = library.code
    if library.adapter:
      self.adapter1 = library.adapter.sequence
    if library.adapter2:
      self.adapter2 = library.adapter2.sequence
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
      fname = os.path.join(lanefile.archive.root_path, library.code, lanefile.filename + ".gz")
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
    self.sample_alias = self.annotations['sample name']
    
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
      sys.exit("Library strategy/type '%s' not supported!\n\n" % (self.library_strategy))

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
      ers = ExternalRecord.objects.get(samples__name=self.sample_alias, repository__name=self.external_repository)
      self.sample_external_record = ers.accession
    except ExternalRecord.DoesNotExist:
      LOGGER.info("No external records for sample \'%s\'.", self.sample_alias)

class EnaXmlObject(object):
  '''A class / container to keep info about generated XML file, its type, alias, validation, submission etc. state'''
  
  def __init__(self, otype, alias, path=None, template_path=None, validated=False, submitted=False):

    if otype not in ['experiment','sample','run']:
      LOGGER.error("Unknown ENA xml submission type \'%s\'. Exiting!", otype)
    
    self.otype = otype # object type. One of following values: experiment, sample, run
    self.alias = alias # experiment, sample or run alias
    self.validated = validated # logical. Shows if XML has been validated
    self.submitted = submitted # logical. Shows if XML has been submitted
    self.accession = None # ENA accession
    fixed_alias = self.fix_alias_for_fname(alias) # alias that has some characters replaced so that it can be used in file names
    
    # check if XML template exists
    self.template = os.path.join(template_path, otype + ".xml")
    if not os.path.exists(self.template):
      LOGGER.error("XML template file %s missig or not accessible!\n\n", self.template)
      sys.exit(1)

    # check if submission XML template exists.
    self.submission_template = os.path.join(template_path, "submission.xml")
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
  
  def __init__(self, lane, study_alias, template_dir="", target_dir="", release_date_str=None, repository_name = 'EBI ENA'):
    '''Takes study alias and optionally, path to xml temlates (template_dir) and targets (dir where xml files should be written).''' 

    self.study_alias = study_alias
    self.template_dir = template_dir
    self.target_dir = target_dir
    self.release_date_str = None
    self.external_repository = repository_name
    
    # Set release date in case not provided, assign it 2 years from now
    if release_date_str is None:
      today = datetime.date.today()
      self.release_date_str = today.replace(year=today.year + 2)
    else:
      self.release_date_str = datetime.strptime(release_date_str , '%Y-%m-%d')
    
    # read annotations from repository
    self.a = EnaAnnotation(lane, self.external_repository, target_dir=self.target_dir)
    self.a.process_annotations()
    self.a.study_alias = study_alias

  def write_experiment_xml(self):
    '''Writes experiment xml file for library. Returns ENA xml info object.'''

    # Build experiment xml info object
    exo = EnaXmlObject('experiment', self.a.experiment_alias, path=self.target_dir, template_path=self.template_dir)

    if self.a.library_external_record is not None:
      exo.submitted = True
      exo.accession = self.a.library_external_record
      return exo
        
    # Check if XML may have already been created
    if os.path.isfile(exo.filename):
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
    exo = EnaXmlObject('run', self.a.run_alias, path=self.target_dir, template_path=self.template_dir)
    
    # Check if XML may have already been created
    if os.path.isfile(exo.filename):
      LOGGER.info("XML for %s \'%s\' already exists: %s", exo.otype, exo.alias, exo.filename)
      # Add filenames to list of datafiles
      for i in range(0,len(self.a.files)):
        exo.datafiles.append(self.a.files[i])
        exo.datamd5s.append(self.a.md5s[i])
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
    exo = EnaXmlObject('sample', self.a.sample_alias, path=self.target_dir, template_path=self.template_dir)

    if self.a.sample_external_record is not None:
      exo.submitted = True
      exo.accession = self.a.sample_external_record
      return exo

    
    # Check if XML may have already been created
    if os.path.isfile(exo.filename):
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
    
    submission = xml.find('SUBMISSION')
    submission.set('alias', exo.alias + "-" + exo.otype + "_submission")        
    submission.set('center_name',self.a.center_name)
    actions = submission.find('ACTIONS')
    i = 1
    for action in actions.findall('ACTION'):
      for child in action:
        if child.tag == 'ADD':
          if child.attrib['schema'] == exo.otype:
            child.set('source', exo.filename)
          else:
            #if child.attrib['schema'] not in ['run','experiment','sample']:
            #  sys.stderr.write('Unsupported ACTION element with schema=\'%s\' in template! Removing the corresponding ACTION element.\n' % child.attrib['schema'])
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


  def upload_xml(self, s_xml, f_xml, s_type, validate=True):
    '''Submits XML to ENA, parses result for failure or of accession. With validation_only being set, submission is directed instead to ENA test server fro XML syntax validation.'''
    # s_xml - submission XML for f_xml file.
    # f_xml - XML file to be submitted.
    # s_type - submission type e.g. ['RUN','ANALYSIS','SAMPLE','EXPERIMENT']
    # validate - if set the submission is directed to ENA test server for XML syntax validation only.
    
    r_fname = s_xml + ".receipt"

    LOGGER.info('Preparing upload of following XMLs to ENA using curl.')
    cmd = "curl -k"
    cmd += " -F \"SUBMISSION=@" + s_xml + "\" -F \"" + s_type.upper() + "=@" + f_xml + "\""
    if validate:
      cmd += " \"https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%20" + self.credentials['username'] + "%20" + self.credentials['password'] + "\""
      r_fname = s_xml + ".test-receipt"
    else:
      cmd += " \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%20" + self.credentials['username'] + "%20" + self.credentials['password'] + "\""
    LOGGER.info(cmd)
    if os.path.isfile(r_fname):
      LOGGER.info("Skipping, already submitted.")
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
      # Parse receipt and exit in case of error
      root = ET.fromstring(stdout)
      submission = root.find('SUBMISSION')
      messages = root.find('MESSAGES')
      errors = messages.findall('ERROR')
      if len(errors) > 0:
        for error in errors:
          LOGGER.error("%s\n", error.text)
          if "already exists as accession" in error.text:
            start = error.text.find('already exists as accession')
            return error.text[(start+27):]
        sys.exit("Invalid XML!\n\n")
      else:
        return submission.attrib['accession']
    else:
      sys.exit("No receipt received for command \"%s\"" % cmd)

  def submit_files(self, file_dict, s_type, test=True):
    '''Submits all files in file dict. File dict contains file names as keys and sublission XML file names for the files as values'''
    for f in file_dict:
      accession = self.upload_xml(file_dict[f], f, s_type, validate=test)
      LOGGER.info("Accession received=%s" % accession)

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

  def __init__(self, study_alias, study_xml, xml_template_dir, credentials_file, release_date_str=None, validate_xml=True, upload_files=False, submit_xml=False, target_dir=""):
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
    
  def submit_lane(self, lane):
    '''Submits data and meta-data for a particular lane_id if not submitted before.'''

    # install XML creator
    ex = EnaXmlCreator(lane, self.study_alias, self.xml_template_dir, self.target_dir, repository_name = self.external_repository)

    # create XMLs
    sample_xml = ex.write_sample_xml()
    experiment_xml = ex.write_experiment_xml()
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
    if not sample_xml.submitted:
      # validate XML
      if self.validate_xml:
        accession = ena_loader.upload_xml(sample_xml.submission_filename, sample_xml.filename, 'SAMPLE', validate=True)
      # submit XML
      if self.submit_xml:
        accession = ena_loader.upload_xml(sample_xml.submission_filename, sample_xml.filename, 'SAMPLE', validate=False)
        erm = ExternalRecordManager(sample_xml.otype, accession = accession, release_date = self.release_date_str, is_public = False, repository_name = self.external_repository)
        # erm.add_sample_obj(sample_xml.alias)
        erm.add_sample_obj_library(ex.a.code)
        erm.add_external_record_to_obj()

    if not experiment_xml.submitted:
      # validate XML
      if self.validate_xml:
        accession = ena_loader.upload_xml(experiment_xml.submission_filename, experiment_xml.filename, 'EXPERIMENT', validate=True)
      # submit XML
      if self.submit_xml:
        accession = ena_loader.upload_xml(experiment_xml.submission_filename, experiment_xml.filename, 'EXPERIMENT', validate=False)
        erm = ExternalRecordManager(experiment_xml.otype, accession = accession, release_date = self.release_date_str, is_public = False, repository_name = self.external_repository)
        erm.add_library_obj(experiment_xml.alias)
        erm.add_external_record_to_obj()
        
    if not run_xml.submitted:
      # validate XML
      if self.validate_xml:
        accession = ena_loader.upload_xml(run_xml.submission_filename, run_xml.filename, 'RUN', validate=True)
      # submit XML
      if self.submit_xml:
        accession = ena_loader.upload_xml(run_xml.submission_filename, run_xml.filename, 'RUN', validate=False)
        erm = ExternalRecordManager(run_xml.otype, accession = accession, release_date = self.release_date_str, is_public = False, repository_name = self.external_repository)
        erm.add_lane_obj(lane.id)
        erm.add_external_record_to_obj()

  def submit_project(self, project):
    '''Submits data and meta-data for all lanes not yet been submitted for a project.'''

    # Find lanes of a project with no association to external repository 'EBI ENA'  
    lanes = Lane.objects.filter(library__projects__name=project).exclude(external_records__repository__name=self.external_repository)
    if lanes.count() == 0:
      LOGGER.error("No unsubmitted lanes for project name=%d!", project)
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

  def submit_lanes_from_file(self, fn):
    '''Submits data and meta-data for a set of lane_ids found in the first column of the file fn'''
    with open(fn, 'rb') as fh:
      for line in fh:
        line = line.rstrip('\n')
        cols = line.split('\t')
        self.submit_lane_with_id(cols[0])
