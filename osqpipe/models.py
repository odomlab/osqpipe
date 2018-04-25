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

'''Core model specification for our repository.'''

import re
import string
import os

from django.db import models
from django.core.exceptions import ValidationError
from django.contrib.auth.models import User

# This stanza to manage REST API authentication tokens (see
# below). Also the project-library ManyToManyField relationship (m2m_changed).
from django.conf import settings
from django.db.models.signals import post_save, m2m_changed
from django.dispatch import receiver
from rest_framework.authtoken.models import Token

import dbarray # https://github.com/ecometrica/django-dbarray
from .managers import ControlledVocabManager, AntibodyManager, FiletypeManager,\
    LibraryManager, LaneManager

from osqutil.config import Config
CONFIG = Config()

####################################################################################
# Data validation functions.

def validate_library_code(value):
  '''Ensures that library codes adopt the basic convention we've had to
  assume to allow correct sorting in the web interface.'''
  value = unicode(value)
  code_re = re.compile('^\w+\d+$')
  if not code_re.match(value):
    raise ValidationError(u'%s is not a valid library code' % value)

####################################################################################
# CV Abstract base class
  
class ControlledVocab(models.Model):

  objects = ControlledVocabManager() # Inherited by subclasses

  _controlled_field = None

  @property
  def controlled_name(self):
    user_key = self.__class__._controlled_field
    if user_key is None:
      raise NotImplementedError()
    else:
      return vars(self)[user_key]

  class Meta:
    abstract = True

####################################################################################
# The concrete model classes themselves.

class Adapter(ControlledVocab):
  code         = models.CharField(max_length=32, unique=True)
  sequence     = models.CharField(max_length=32, null=True, blank=True)
  protocol     = models.CharField(max_length=32)

  _controlled_field = 'code'
  
  def __unicode__(self):
    return self.code

  class Meta:
    db_table = u'adapter'
    ordering = ['code']

class Program(models.Model):
  program      = models.CharField(max_length=128)
  version      = models.CharField(max_length=128)
  options      = models.CharField(max_length=256, null=True, blank=True)
  files        = models.CharField(max_length=256, null=True, blank=True)
  type         = models.CharField(max_length=128) # FIXME we need another CV table here
  description  = models.TextField(null=True, blank=True)
  date         = models.DateField(auto_now_add=True)
  current      = models.BooleanField(default=True)

  def __unicode__(self):
    return "%s (version %s)" % (self.program, self.version)

  class Meta:
    db_table = u'program'
    unique_together = ('program', 'version')
    ordering = ['program', 'version']

class Machine(ControlledVocab):
  code         = models.CharField(max_length=32, unique=True)
  platform     = models.CharField(max_length=32)
  name         = models.CharField(max_length=32)
  
  _controlled_field = 'code'
  
  def __unicode__(self):
    return self.code

  class Meta:
    db_table = u'machine'
    verbose_name_plural = 'machines'
    ordering = ['code']

class Species(models.Model):
  scientific_name = models.CharField(max_length=255, db_column='sciname', unique=True)
  common_name     = models.CharField(max_length=255, db_column='commonname', null=True, blank=True)
  accession       = models.CharField(max_length=32, unique=True)

  def __unicode__(self):
    return self.scientific_name

  class Meta:
    db_table = u'species'
    verbose_name_plural = 'species'
    ordering = ['scientific_name']

class Genome(ControlledVocab):
  code         = models.CharField(max_length=32, unique=True)
  blastdb      = models.CharField(max_length=255, null=True, blank=True)
  fasta        = models.CharField(max_length=255, null=True, blank=True)
  fasta_md5sum = models.CharField(max_length=32, null=True, blank=True)
  notes        = models.TextField(null=True, blank=True)
  version      = models.CharField(max_length=255, null=True, blank=True)
  url          = models.CharField(max_length=256, null=True, blank=True)
  species      = models.ForeignKey(Species, on_delete=models.PROTECT)

  _controlled_field = 'code'
  
  @property
  def fasta_path(self):
    '''
    Returns the expected path to the fasta file for a given
    genome. This is distinct from the fasta files associated with
    various indices (bwa, tophat, etc.) which reside in their own
    subdirectories.
    '''
    sciname = str(self.species.scientific_name)
    sciname = sciname.replace(" ", "_")
    sciname = sciname.lower()
    fasta   = os.path.join(CONFIG.clustergenomedir, sciname,
                           self.code, "%s.fa" % self.code)

    return fasta

  def __unicode__(self):
    return self.code

  class Meta:
    db_table = u'genome'
    ordering = ['code']

class Restrictome(ControlledVocab):
  enzyme       = models.CharField(max_length=128)
  sequence     = models.CharField(max_length=128, null=True, blank=True)
  filename     = models.CharField(max_length=1024, unique=True)
  date         = models.DateField(auto_now_add=True)
  program      = models.ForeignKey(Program, on_delete=models.PROTECT, help_text="Program that was used for creating restrictome.")
  genome       = models.ForeignKey(Genome, on_delete=models.PROTECT,
                                   help_text="The genome against which the restrictome was generated.")
  def __unicode__(self):
    return "%s %s" % (self.enzyme, self.genome.code)

  class Meta:
    db_table = u'restrictome'
    ordering = ['enzyme','genome']
   
class Tissue(ControlledVocab):
  name         = models.CharField(max_length=255, unique=True)
  description  = models.TextField(null=True, blank=True)

  _controlled_field = 'name'

  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'tissue'
    ordering = ['name']

class Antibody(ControlledVocab):
  name         = models.CharField(max_length=255)
  lot_number   = models.CharField(max_length=64, default='unknown')
  description  = models.TextField(null=True, blank=True)

  objects      = AntibodyManager()

  _controlled_field = 'name'
  
  def __unicode__(self):
    return "%s (lot: %s)" % (self.name, self.lot_number)

  class Meta:
    db_table = u'antibody'
    unique_together = ('name', 'lot_number')
    verbose_name_plural = 'antibodies'
    ordering = ['name', 'lot_number']
    
class Factor(ControlledVocab):
  name         = models.CharField(max_length=255, unique=True)
  description  = models.TextField(null=True, blank=True)

  _controlled_field = 'name'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'factor'
    ordering = ['name']

class Strain(ControlledVocab):
  name         = models.CharField(max_length=255, unique=True)
  description  = models.TextField(null=True, blank=True)

  _controlled_field = 'name'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'strain'
    ordering = ['name']

class Sex(ControlledVocab):
  name         = models.CharField(max_length=32, unique=True)

  _controlled_field = 'name'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'sex'
    ordering = ['name']

class Libtype(ControlledVocab):
  code         = models.CharField(max_length=32, unique=True)
  name         = models.CharField(max_length=32, unique=True)
  description  = models.CharField(max_length=255, null=True, blank=True)

  _controlled_field = 'code'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'libtype'
    ordering = ['code']

class Condition(ControlledVocab):
  name         = models.CharField(max_length=255, unique=True)
  description  = models.TextField(null=True, blank=True)

  _controlled_field = 'name'

  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'condition'
    ordering = ['name']

class Linkerset(ControlledVocab):
  name         = models.CharField(max_length=255, unique=True)
  fivep        = models.CharField(max_length=255)
  threep       = models.CharField(max_length=255)

  _controlled_field = 'name'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'linkerset'
    ordering = ['name']

class Project(ControlledVocab):
  code         = models.CharField(max_length=32, unique=True)
  name         = models.CharField(max_length=32, unique=True)
  description  = models.TextField(null=True, blank=True)
  shortnames   = models.BooleanField(default=True)
  filtered     = models.BooleanField(default=True)
  lab          = models.CharField(max_length=32)
  people       = models.ManyToManyField(User, db_table='project_users')
  is_frozen    = models.BooleanField(default=False)
  # libraries relationship dealt with below.

  _controlled_field = 'code'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'project'
    ordering = ['name']

class Source(models.Model):
  '''
  The original source organism from which tissue was taken.
  '''
  name         = models.CharField(max_length=128, unique=True)
  strain       = models.ForeignKey(Strain, on_delete=models.PROTECT, null=True, blank=True)
  sex          = models.ForeignKey(Sex, on_delete=models.PROTECT, null=True, blank=True)
  date_of_birth = models.DateField(null=True, blank=True)
  date_of_death = models.DateField(null=True, blank=True)
  mother       = models.ForeignKey('Source', on_delete=models.PROTECT, null=True, blank=True, related_name='child_as_mother')
  father       = models.ForeignKey('Source', on_delete=models.PROTECT, null=True, blank=True, related_name='child_as_father')
  species      = models.ForeignKey(Species, on_delete=models.PROTECT, null=True, blank=True)
  comment      = models.TextField(null=True, blank=True)

  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'source'
    ordering = ['name']

class DoseUnit(ControlledVocab):
  '''
  Chemical concentration, radiation or other unit of dosage.
  '''
  name         = models.CharField(max_length=32, unique=True)
  description  = models.CharField(max_length=128)

  _controlled_field = 'name'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'dose_unit'
    ordering = ['name']

class TreatmentAgent(ControlledVocab):
  '''
  Chemical or other agent used to treat a Source organism prior to
  taking a Sample.
  '''
  name         = models.CharField(max_length=64, unique=True)
  description  = models.CharField(max_length=256)
  accession    = models.CharField(max_length=32)

  _controlled_field = 'name'

  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'treatment_agent'
    ordering = ['name']

class SourceTreatment(models.Model):
  '''
  Each Source (i.e., originating organism) may be treated multiple
  times in multiple ways. Note that dose is deliberately coded as
  CharField to allow for maximum flexibility.
  '''
  source       = models.ForeignKey(Source, on_delete=models.PROTECT)
  date         = models.DateField()
  agent        = models.ForeignKey(TreatmentAgent, on_delete=models.PROTECT)
  dose         = models.CharField(max_length=128, null=True, blank=True)
  dose_unit    = models.ForeignKey(DoseUnit, on_delete=models.PROTECT,
                                   null=True, blank=True)

  def __unicode__(self):
    retval = "%s %s" % (self.date, self.agent)
    if self.dose is not None:
      retval += " (%s" % self.dose
      if self.dose_unit is not None:
        retval += " %s)" % self.dose_unit
      else:
        retval += ")"
    return retval

  class Meta:
    db_table = u'source_treatment'
    ordering = ['date', 'agent']
    unique_together = ['source', 'date', 'agent']

# Technically this should perhaps be a ControlledVocab subclass, but
# at the moment CV only uses the _controlled_name field to search,
# without reference to category. This is therefore something a bit
# different. Re-examine this if we ever need fuzzy search matching.
class Characteristic(models.Model):
  '''
  A class designed to store category-value terms. These should be
  characteristics with discrete, defined values (i.e. categorical
  variables, *not* continuous variables which might be better modelled
  as a Measurement class.
  '''
  category     = models.CharField(max_length=32)
  value        = models.CharField(max_length=32)

  def __unicode__(self):
    return "%s: %s" % (self.category, self.value)

  class Meta:
    db_table = u'characteristic'
    ordering = ['category','value']
    unique_together = ['category','value']

class SizeUnit(ControlledVocab):
  '''
  A unit of size. This is deliberately vague as it may be a length,
  volume, or conceivably even a mass unit.
  '''
  name         = models.CharField(max_length=32, unique=True)
  description  = models.CharField(max_length=128)

  _controlled_field = 'name'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'size_unit'
    ordering = ['name']

class ExternalRepository(ControlledVocab):
  name         = models.CharField(max_length=32, unique=True)

  _controlled_field = 'name'

  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'external_repository'
    ordering = ['name']

class ExternalRecord(ControlledVocab):
  accession    = models.CharField(max_length=32, unique=True)
  repository   = models.ForeignKey(ExternalRepository, on_delete=models.PROTECT)
  is_public    = models.BooleanField(default=False)
  release_date = models.DateField()

  _controlled_field = 'accession'
  
  def __unicode__(self):
    return "%s: %s" % (self.repository, self.accession)

  class Meta:
    db_table = u'external_record'
    ordering = ['accession']

    
class Sample(models.Model):
  '''
  A tissue sample taken from a Source organism.
  '''
  name         = models.CharField(max_length=128)
  tissue       = models.ForeignKey(Tissue, on_delete=models.PROTECT)
  source       = models.ForeignKey(Source, on_delete=models.PROTECT)
  characteristics = models.ManyToManyField(Characteristic, db_table='sample_characteristic', related_name='samples')
  size         = models.DecimalField(max_digits=5, decimal_places=2, null=True, blank=True)
  size_unit    = models.ForeignKey(SizeUnit, on_delete=models.PROTECT,
                                   null=True, blank=True)
  comment      = models.TextField(null=True, blank=True)
  external_records = models.ManyToManyField(ExternalRecord, db_table='sample_external_record', related_name='samples')

  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'sample'
    ordering = ['name']
    unique_together = ['name', 'tissue']

class Library(models.Model):
  code         = models.CharField(max_length=128, unique=True,
                                  validators=[validate_library_code])
  genome       = models.ForeignKey(Genome, on_delete=models.PROTECT,
                                   help_text="The genome against which the sequence"
                                   + " data from this library should be aligned.")
  sample       = models.ForeignKey(Sample, on_delete=models.PROTECT)
  antibody     = models.ForeignKey(Antibody, on_delete=models.PROTECT, null=True, blank=True)
  factor       = models.ForeignKey(Factor, on_delete=models.PROTECT, null=True, blank=True)
  condition    = models.ForeignKey(Condition, on_delete=models.PROTECT, null=True, blank=True)
  bad          = models.BooleanField(default=False)
  release_worthy = models.BooleanField(default=False)
  projects     = models.ManyToManyField(Project, db_table='library_project', related_name='libraries')
  libtype      = models.ForeignKey(Libtype, on_delete=models.PROTECT)
  barcode      = models.CharField(max_length=32, null=True, blank=True)
  linkerset    = models.ForeignKey(Linkerset, on_delete=models.PROTECT, null=True, blank=True)
  chipsample   = models.CharField(max_length=255, null=True, blank=True)
  paired       = models.BooleanField(default=False)
  adapter      = models.ForeignKey(Adapter, on_delete=models.PROTECT, null=True, blank=True, related_name='libraries')
  adapter2     = models.ForeignKey(Adapter, on_delete=models.PROTECT, null=True, blank=True, related_name='libraries2')
  platecode    = models.CharField(max_length=32, null=True, blank=True)
  platecol     = models.IntegerField(null=True, blank=True)
  platerow     = models.CharField(max_length=1,
                                  choices=[('A','A'),('B','B'),('C','C'),('D','D'),
                                           ('E','E'),('F','F'),('G','G'),('H','H')],
                                  null=True, blank=True)
  comment      = models.TextField(null=True, blank=True)
  external_records = models.ManyToManyField(ExternalRecord, db_table='library_external_record', related_name='libraries')
  objects      = LibraryManager()

  @property
  def filename_tag(self):
    """
    On-the-fly replacement for the old library description field, with
    the library code included for completeness.
    """
    fac = self.factor.name   if self.factor     else 'unk'
    ant = self.antibody.name if self.antibody   else 'unk'
    ind = self.sample.name   if self.sample.name != self.code else ''
    sta = self.sample.source.strain.name if self.sample.source.strain  else ''

    tis = self.sample.tissue.name
    if self.condition is not None:
      tis = "%s_%s" % (tis, self.condition)
    gen = self.genome.code

    # Here we pull out unnecessary duplication of the tissue name
    # which may have been introduced when generating a new sample name
    # (see osqpipe.pipeline.library).
    ind = string.replace(ind, tis, '')

    # Try to strip out not just trailing space, but some other characters we might use at a later date.
    ind = re.sub('[- _/;:.,]$', '', ind)

    return("%s_%s_%s_%s_%s%s%s"
           % (self.code, fac, tis, ant, gen, sta, ind))

  def __unicode__(self):
    return self.code

  class Meta:
    db_table = u'library'
    verbose_name_plural = 'libraries'
    ordering = ['code']
    
class Facility(ControlledVocab):
  code         = models.CharField(max_length=10, unique=True)
  name         = models.CharField(max_length=255, unique=True)
  description  = models.TextField(null=True, blank=True)

  _controlled_field = 'code'
  
  def __unicode__(self):
    return self.code

  class Meta:
    db_table = u'facility'
    verbose_name_plural = 'facilities'
    ordering = ['code']

class Status(ControlledVocab):
  code         = models.CharField(max_length=32, unique=True)
  description  = models.CharField(max_length=128)
  colour       = models.CharField(max_length=32, default='#FFFFFF')
  lanerelevant    = models.BooleanField(default=False)
  libraryrelevant = models.BooleanField(default=False)
  sortcode     = models.IntegerField(default=0)
  authority    = models.ForeignKey(Facility, on_delete=models.PROTECT, null=True, blank=True)

  _controlled_field = 'code'
  
  def __unicode__(self):
    if self.authority is not None:
      return "%s: %s" % (self.authority.name, self.code)
    else:
      return self.code

  @property
  def display_tag(self):
    return unicode(self)
    
  class Meta:
    db_table = u'status'
    verbose_name_plural = 'status flags'
    ordering = ['code']

class Lane(models.Model):
  library      = models.ForeignKey(Library, on_delete=models.PROTECT)
  machine      = models.ForeignKey(Machine, on_delete=models.PROTECT)
  flowcell     = models.CharField(max_length=32)
  rundate      = models.DateField()
  reads        = models.IntegerField(null=True, blank=True)
  passedpf     = models.IntegerField(null=True, blank=True)
  lanenum      = models.IntegerField()
  flowlane     = models.IntegerField()
  paired       = models.BooleanField(default=False)
  readlength   = models.IntegerField(null=True, blank=True)
  mapped       = models.IntegerField(null=True, blank=True)
  facility     = models.ForeignKey(Facility, on_delete=models.PROTECT)
  seqsamplepf  = models.TextField(blank=True)
  seqsamplebad = models.TextField(blank=True)
  qualmeanpf   = dbarray.FloatArrayField(null=True)
  qualstdevpf  = dbarray.FloatArrayField(null=True)
  qualmean     = dbarray.FloatArrayField(null=True)
  qualstdev    = dbarray.FloatArrayField(null=True)
  summaryurl   = models.CharField(max_length=1024, null=True, blank=True)
  genomicssampleid = models.CharField(max_length=32, null=True, blank=True)
  usersampleid     = models.CharField(max_length=1024, null=True, blank=True)
  notes        = models.TextField(null=True, blank=True)
  failed       = models.BooleanField(default=False)
  runnumber    = models.CharField(null=True, blank=True, max_length=255)
  status       = models.ForeignKey(Status, on_delete=models.PROTECT)
  external_records = models.ManyToManyField(ExternalRecord, db_table='lane_external_record', related_name='lanes')

  objects      = LaneManager()

  @property
  def total_passedpf(self):
    if self.paired and self.passedpf is not None:
      return self.passedpf * 2
    else:
      return self.passedpf

  @property
  def total_reads(self):
    if self.paired and self.reads is not None:
      return self.reads * 2
    else:
      return self.reads

  @property
  def name(self):
    return "%s%02d" % (self.facility, self.lanenum)

  @property
  def notes_formatted(self):
    if self.notes is None:
      return ''
    else:
      return "\n".join(self.notes.split(";"))

  @property
  def model_parent(self):
    return self.library

  @property
  def is_published(self):
    return any([ x.is_public for x in self.external_records.all() ])

  @property
  def public_records(self):
    return [ x for x in self.external_records.all() if x.is_public ]

  def __unicode__(self):
    return "%s_%s%02d" % (self.library, self.facility, self.lanenum)

  class Meta:
    db_table = u'lane'
    unique_together = ('library', 'lanenum', 'facility')
    ordering = ['library']

####################################################################################
# Data processing abstract base class.

class DataProcess(models.Model):

  # Links to DataProvenance. This has to be a concrete base class
  # (i.e. having its own database table) to allow that link to be
  # created.

  def __unicode__(self):
    return ", ".join([str(x) for x in self.provenance.all().order_by('rank_index')])

  class Meta:
    db_table = u'data_process'

class DataProvenance(models.Model):

  rank_index   = models.IntegerField()
  program      = models.ForeignKey(Program, on_delete=models.PROTECT)
  parameters   = models.CharField(max_length=255)
  data_process = models.ForeignKey(DataProcess, related_name='provenance')

  def __unicode__(self):
    if len(self.parameters) > 0:
      return "%s %s" % (self.program, self.parameters)
    else:
      return unicode(self.program)

  class Meta:
    db_table = u'data_provenance'
    unique_together = ('data_process', 'rank_index')

class Alignment(DataProcess):
  genome       = models.ForeignKey(Genome, on_delete=models.PROTECT)
  total_reads  = models.IntegerField()
  mapped       = models.IntegerField()
  munique      = models.IntegerField(null=True, blank=True)
  headtrim     = models.IntegerField(default=0, blank=True)
  tailtrim     = models.IntegerField(default=0, blank=True)
  lane         = models.ForeignKey(Lane, on_delete=models.PROTECT)

  @property
  def mapped_percent(self):
    if self.lane.total_passedpf in (0,None):
      return None
    else:
      return round(100*(float(self.mapped)/self.lane.total_passedpf), 1)

  @property
  def munique_percent(self):
    if self.lane.total_passedpf in (0,None):
      return None
    else:
      return round(100*(float(self.munique)/self.lane.total_passedpf), 1)

  def __unicode__(self):
    provenance = ", ".join([str(x) for x in self.provenance.all().order_by('rank_index')])
    return "%s : %s (%s)" % (self.lane, self.genome, provenance)

  class Meta:
    db_table = u'alignment'

class MergedAlignment(DataProcess):
  alignments = models.ManyToManyField(Alignment)

  @property
  def genome(self):
    self.full_clean() # Ensures only one genome linked via alignments.
    # FIXME use the first() method once we've migrated to django >= 1.6
    return self.alignments.all()[:1].get().genome

  # Custom model validation here.
  def clean(self):
    if self.alignments.count() == 0:
      raise ValidationError(\
        {'alignments' :
           'MergedAlignment is linked to zero source Alignments.'})
    genomes = self.alignments.values_list('genome', flat=True).distinct()
    if genomes.count() != 1:
      raise ValidationError(\
        {'alignments' :
           'MergedAlignment links Alignments against multiple different genome builds.'})

  def __unicode__(self):
    return "%s (%s)" % (";".join([str(x.lane) for x in self.alignments.all()]), self.genome)

  class Meta:
    db_table = u'merged_alignment'

class LaneQC(DataProcess):

  lane        = models.ForeignKey(Lane, on_delete=models.PROTECT)

  def __unicode__(self):
    provenance = ", ".join([str(x) for x in self.provenance.all().order_by('rank_index')])
    return "%s %s" % (self.lane, provenance)

  class Meta:
    db_table = u'lane_qc'
    verbose_name = u'Lane QC'

class QCValue(models.Model):
  name         = models.CharField(max_length=32)
  value        = models.CharField(max_length=32)
  laneqc       = models.ForeignKey(LaneQC, on_delete=models.CASCADE, related_name='qc_values')

  def __unicode__(self):
    return "%s %s=%s" % (self.laneqc, self.name, self.value)

  class Meta:
    db_table        = u'qc_value'
    verbose_name    = u'QC Value'
    unique_together = ('laneqc', 'name')

class AlignmentQC(DataProcess):

  alignment        = models.ForeignKey(Alignment, on_delete=models.PROTECT)

  def __unicode__(self):
    provenance = ", ".join([str(x) for x in self.provenance.all().order_by('rank_index')])
    return "%s %s" % (self.alignment, provenance)

  class Meta:
    db_table = u'alignment_qc'
    verbose_name = u'Alignment QC'

class ArchiveLocation(models.Model):
  name         = models.CharField(max_length=32, unique=True)
  root_path    = models.CharField(max_length=1024, unique=True)
  host         = models.CharField(max_length=1024, unique=False)
  host_port    = models.CharField(max_length=8, unique=False)
  host_path    = models.CharField(max_length=1024, unique=False)
  host_user    = models.CharField(max_length=128, unique=False)
  host_delete_timelag = models.IntegerField(editable=True, null=True)

  def __unicode__(self):
    return "%s (%s)" % (self.name, self.root_path)

  class Meta:
    db_table        = u'archive_location'

class Filetype(ControlledVocab):
  code         = models.CharField(max_length=10, unique=True)
  name         = models.CharField(max_length=32, unique=True)
  description  = models.TextField(null=True, blank=True)
  suffix       = models.CharField(max_length=32, blank=True, unique=True)
  gzip         = models.BooleanField(default=True)

  objects      = FiletypeManager()

  _controlled_field = 'code'
  
  def __unicode__(self):
    return self.name

  class Meta:
    db_table = u'filetype'
    ordering = ['code']

####################################################################################
# Datafile abstract base class.

class Datafile(models.Model):
  filename     = models.CharField(max_length=1024, unique=True)
  checksum     = models.CharField(max_length=128)
  filetype     = models.ForeignKey(Filetype, on_delete=models.PROTECT)
  description  = models.TextField(null=True, blank=True)
  date         = models.DateField(auto_now_add=True)
  archive      = models.ForeignKey(ArchiveLocation, on_delete=models.PROTECT,
                                    null=True, blank=True)
  archive_date = models.DateField(auto_now_add=False, null=True, blank=True)

  @property
  def libcode(self):
    raise NotImplementedError();

  @property
  def repository_root(self):
    if self.archive is None:
      return CONFIG.repositorydir
    else:
      return self.archive.root_path

  @property
  def filename_on_disk(self):
    fname   = self.filename
    if self.filetype.gzip:
      fname += CONFIG.gzsuffix
    return fname

  @property
  def repository_file_path(self):
    '''
    The current path to the file on disk.
    '''
    return os.path.join(self.repository_root, self.libcode, self.filename_on_disk)

  @property
  def original_repository_file_path(self):
    '''
    The original repository path to the file on disk, prior to any
    archival.
    '''
    return os.path.join(CONFIG.repositorydir, self.libcode, self.filename_on_disk)
        
  def __unicode__(self):
    return self.filename

  class Meta:
    abstract = True

class Lanefile(Datafile):
  lane         = models.ForeignKey(Lane, on_delete=models.PROTECT)
  pipeline     = models.CharField(max_length=128, default='chipseq') # FIXME

  @property
  def libcode(self):
    return self.lane.library.code

  class Meta:
    db_table = u'lanefile'
    ordering = ['filename']

class QCfile(Datafile):
  laneqc      = models.ForeignKey(LaneQC, on_delete=models.PROTECT)

  @property
  def libcode(self):
    return self.laneqc.lane.library.code

  class Meta:
    db_table = u'qcfile'
    verbose_name = u'QC file'
    ordering = ['filename']

class Alnfile(Datafile):
  alignment    = models.ForeignKey(Alignment, on_delete=models.PROTECT)

  @property
  def libcode(self):
    return self.alignment.lane.library.code

  class Meta:
    db_table = u'alnfile'
    ordering = ['filename']

class AlnQCfile(Datafile):
  alignmentqc      = models.ForeignKey(AlignmentQC, on_delete=models.PROTECT)

  @property
  def libcode(self):
    return self.alignmentqc.alignment.lane.library.code

  class Meta:
    db_table = u'alnqcfile'
    verbose_name = u'Alignment QC file'
    ordering = ['filename']

class MergedAlnfile(Datafile):
  alignment    = models.ForeignKey(MergedAlignment, on_delete=models.PROTECT)

  @property
  def repository_file_path(self):
    fname   = self.filename
    if self.filetype.gzip:
      fname += CONFIG.gzsuffix

    return os.path.join(self.repository_root, CONFIG.merged_alignments_dir, fname)
    
  class Meta:
    db_table = u'merged_alnfile'
    ordering = ['filename']

class HistologyImagefile(Datafile):
  sample    = models.ForeignKey(Sample, on_delete=models.PROTECT)
  block     = models.CharField(max_length=32)
  batch     = models.CharField(max_length=32)

  @property
  def repository_file_path(self):
    fname   = self.filename
    if self.filetype.gzip:
      fname += CONFIG.gzsuffix

    return os.path.join(self.repository_root, CONFIG.histology_image_dir, fname)
    
  class Meta:
    db_table = u'histology_imagefile'
    ordering = ['filename']

# FIXME could also be deprecated, might conceivably be useful though.
class LibraryNameMap(models.Model):
  libname      = models.CharField(max_length=128)
  limsname     = models.CharField(max_length=128, unique=True)

  def __unicode__(self):
    return "%s == %s" % (self.limsname, self.libname)

  class Meta:
    db_table = u'library_name_map'

class Peakcalling(DataProcess):
  code         = models.CharField(max_length=128, unique=True)
  factor_align = models.ForeignKey(Alignment, on_delete=models.PROTECT,
                                   related_name='+', db_column='flane_id')
  input_align  = models.ForeignKey(Alignment, on_delete=models.PROTECT,
                                   related_name='+', db_column='ilane_id')
  description  = models.TextField(null=True, blank=True)
  date         = models.DateField(auto_now_add=True)

  def __unicode__(self):
    provenance = ", ".join([str(x) for x in self.provenance.all().order_by('rank_index')])
    return "%s: %s (input: %s) using %s" % (self.code, self.factor_align,
                                            self.input_align, provenance)

  class Meta:
    db_table = u'peakcalls'

class Peakfile(Datafile):
  peakcalling  = models.ForeignKey(Peakcalling, on_delete=models.PROTECT,
                                   db_column='peakcalls_id')

  @property
  def libcode(self):
    # We need to decide where we're putting these FIXME.
    raise NotImplementedError()
  
  class Meta:
    db_table = u'peakfile'
    ordering = ['filename']

####################################################################################
# An unmanaged class to allow us to access the library_extra table
# which splits library names into prefix+suffix for sorting purposes.
    
class LibraryExtra(models.Model):
  library = models.OneToOneField(Library, related_name='extra',
                                 editable=False, on_delete=models.DO_NOTHING)
  code_text_prefix = models.CharField(max_length=128, editable=False)
  code_numeric_suffix = models.IntegerField(editable=False, null=True)

  class Meta:
    db_table = u'library_extra'
    managed  = False

####################################################################################
# Catch the User objects' post-save signals to create a Token for the REST API.
    
@receiver(post_save, sender=settings.AUTH_USER_MODEL)
def create_auth_token(sender, instance=None, created=False, **kwargs):
  if created:
    Token.objects.create(user=instance)

####################################################################################
# Manage the project-library relationship for projects which are data-frozen.

@receiver(m2m_changed, sender=Library.projects.through)
def protect_frozen_m2m(sender, instance, action, reverse, model, pk_set, **kwargs):
  '''
  Signal receiver to monitor data-frozen projects and prevent
  addition or removal of libraries.
  '''
  # We trigger exceptions before addition/removal/clearance, not after.
  if action in ('pre_add', 'pre_remove', 'pre_clear'):

    # instance can be either Project or Library.
    if reverse:

      # instance is project; model is library.
      if instance.is_frozen:
        raise ValidationError("Attempt to change a frozen project instance.")

    elif pk_set is not None:

      # instance is library; model is project.
      for pkid in list(pk_set):
        if model.objects.get(id=pkid).is_frozen:
          raise ValidationError("Attempt to change a frozen project instance.")
