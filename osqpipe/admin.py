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

'''Configuration for the auto-generated django admin web interface.'''

from .models import *
from django.contrib import admin
from django.core.urlresolvers import reverse
from django import forms

# Admin-level display config is held in each *Admin class; we register
# all model components except for LibraryNameIgnore which I suspect
# could be deprecated.

#############################################
@admin.register(Adapter)
class AdapterAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'sequence', 'protocol')
  search_fields = ('code', 'sequence', 'protocol')

#############################################
#FIXME we need to handle provenance (program/parameters) here
@admin.register(Alignment)
class AlignmentAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'lane_link', 'genome')
  search_fields = ('lane__library__code', 'genome__code')
  readonly_fields = ('lane',)
  fields = ('lane', 'genome',
            'mapped', 'munique', 'total_reads', 'headtrim', 'tailtrim')

  def lane_link(self, obj):
    url = reverse('admin:osqpipe_lane_change', args=(obj.lane.pk,))
    return '<a href="%s">%s</a>' % (url, obj.lane)
  lane_link.allow_tags = True

#############################################
@admin.register(AlignmentQC)
class AlignmentQCAdmin(admin.ModelAdmin):
  # FIXME provenance here also
  list_display = ('__unicode__', 'alignment_link')
  search_fields = ('alignment__lane__library__code',)
  readonly_fields = ('alignment',)
  fields = ('alignment',)

  def alignment_link(self, obj):
    url = reverse('admin:osqpipe_alignment_change', args=(obj.alignment.pk,))
    return '<a href="%s">%s</a>' % (url, obj.alignment)
  alignment_link.allow_tags = True

#############################################
@admin.register(Alnfile)
class AlnfileAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'alignment_link', 'filetype')

  search_fields   = ('filename', 'filetype__name')
  readonly_fields = ('alignment', 'date')
  
  fields = ('filename', 'checksum', 'filetype',
            'alignment', 'date', 'description')

  def alignment_link(self, obj):
    url = reverse('admin:osqpipe_alignment_change', args=(obj.alignment.pk,))
    return '<a href="%s">%s</a>' % (url, obj.alignment)
  alignment_link.allow_tags = True

#############################################
@admin.register(Antibody)
class AntibodyAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'lot_number', 'description')
  search_fields = ('name', 'lot_number')

#############################################
@admin.register(Facility)
class FacilityAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'name', 'description')

#############################################
@admin.register(Factor)
class FactorAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name',)

#############################################
@admin.register(Filetype)
class FiletypeAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'name', 'suffix', 'description')

#############################################
@admin.register(Genome)
class GenomeAdmin(admin.ModelAdmin):
  list_display  = ('__unicode__', 'version')
  
  search_fields = ('code', 'species__common_name', 'species__scientific_name')

  fields = ('code', 'species', 'fasta',
            'fasta_md5sum', 'url', 'notes', 'version', 'blastdb')

#############################################
@admin.register(Restrictome)
class RestrictomeAdmin(admin.ModelAdmin):
  list_display  = ('enzyme', 'genome')
  
  search_fields = ('enzyme', 'genome__code', 'sequence')

  fields = ('genome', 'enzyme', 'sequence', 'program', 'filename')

#############################################
@admin.register(Species)
class SpeciesAdmin(admin.ModelAdmin):
  list_display  = ('scientific_name', 'common_name', 'accession')
  
  search_fields = ('common_name', 'scientific_name')

  fields = ('scientific_name', 'common_name', 'accession')

#############################################
@admin.register(Lane)
class LaneAdmin(admin.ModelAdmin):
  list_display  = ('library', 'facility', 'lanenum',
                   'flowcell', 'runnumber', 'rundate')

  search_fields = ('library__code', 'facility__name', 'flowcell', 'runnumber')

  readonly_fields = ('external_records',)

  fieldsets = [
    (None,
     {'fields': ['library', 'facility', 'lanenum']}),
    ('Run Data',
     {'fields': ['flowcell', 'flowlane', 'runnumber', 'machine', 'rundate', 'paired', 'readlength']}),
    ('LIMS Record',
     {'fields': ['genomicssampleid', 'usersampleid']}),
    ('QC Data',
     {'fields': ['reads', 'passedpf', 'mapped', 'failed', 'status', 'summaryurl']}),
    ('Other Info',
     {'fields': ['external_records', 'notes'] })
    ]

#############################################
@admin.register(LaneQC)
class LaneQCAdmin(admin.ModelAdmin):
  # FIXME provenance here also
  list_display = ('__unicode__', 'lane_link')
  search_fields = ('lane__library__code',)
  readonly_fields = ('lane',)
  fields = ('lane',)

  def lane_link(self, obj):
    url = reverse('admin:osqpipe_lane_change', args=(obj.lane.pk,))
    return '<a href="%s">%s</a>' % (url, obj.lane)
  lane_link.allow_tags = True

#############################################
@admin.register(Lanefile)
class LanefileAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'lane_link', 'filetype')

  search_fields   = ('filename', 'filetype__name')
  readonly_fields = ('lane', 'pipeline')

  fields = ('filename', 'checksum', 'filetype',
            'lane', 'pipeline', 'description')
    
  def lane_link(self, obj):
    url = reverse('admin:osqpipe_lane_change', args=(obj.lane.pk,))
    return '<a href="%s">%s</a>' % (url, obj.lane)
  lane_link.allow_tags = True

#############################################
@admin.register(Characteristic)
class CharacteristicAdmin(admin.ModelAdmin):
  list_display  = ('category', 'value')
  
  search_fields = ('category','value')

  fields = ('category','value')

#############################################
@admin.register(Sample)
class SampleAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'source', 'tissue')

  search_fields  = ('name', 'source__name', 'tissue__name', 'comment')

  readonly_fields = ('source',)
  fields = ('name', 'source', 'tissue', 'comment')

#############################################
@admin.register(Source)
class SourceAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'strain', 'sex')

  search_fields  = ('name', 'strain__name', 'sex__name', 'species__scientific_name')

  fields = ('name', 'species', 'strain', 'sex', 'date_of_birth',
            'date_of_death', 'mother', 'father', 'comment')

#############################################
@admin.register(SourceTreatment)
class SourceTreatmentAdmin(admin.ModelAdmin):
  list_display       = ('source', 'agent', 'dose', 'dose_unit')

  search_fields  = ('source', 'agent')

  fields = ('agent', 'date', 'dose', 'dose_unit')

#############################################
class LibraryAdminForm(forms.ModelForm):

  '''A custom form allowing us to override the default height of
  multi-select boxes.'''

  class Meta:
    widgets = {
      'projects': forms.SelectMultiple(attrs={ 'size': 15 })
      }

@admin.register(Library)
class LibraryAdmin(admin.ModelAdmin):
  form = LibraryAdminForm
  
  list_display       = ('__unicode__', 'genome_link', 'libtype', 'chipsample',
                        'sample_link', 'factor', 'antibody')

  search_fields  = ('code', 'genome__code', 'libtype__name', 'sample__source__sex__name',
                    'sample__source__strain__name', 'sample__name', 'chipsample',
                    'factor__name', 'sample__tissue__name', 'antibody__name')

  fieldsets = [
    (None,
     {'fields': ['code'] }),
    ('Sample Annotation',
     {'fields': ['sample'] }),
    ('Library Details',
     {'fields': ['factor', 'antibody'] }),
    ('Technical Info',
     {'fields': ['libtype', 'genome', 'paired',
                 'adapter', 'adapter2', 'linkerset', 'barcode', 'bad']}),
    ('Project Info',
     {'fields': ['projects', 'chipsample', 'comment']})
    ]

  def sample_link(self, obj):
    url = reverse('admin:osqpipe_sample_change', args=(obj.sample.pk,))
    return '<a href="%s">%s</a>' % (url, obj.sample)
  sample_link.allow_tags = True
  def genome_link(self, obj):
    url = reverse('admin:osqpipe_genome_change', args=(obj.genome.pk,))
    return '<a href="%s">%s</a>' % (url, obj.genome)
  genome_link.allow_tags = True

#############################################
@admin.register(LibraryNameMap)
class LibraryNameMapAdmin(admin.ModelAdmin):
  list_display  = ('limsname', 'libname')
  fields        = ('limsname', 'libname')
  search_fields = ('limsname', 'libname')

#############################################
@admin.register(Libtype)
class LibtypeAdmin(admin.ModelAdmin):
  list_display = ('name', 'code', 'description')
  fields       = ('name', 'code', 'description')

#############################################
@admin.register(Linkerset)
class LinkersetAdmin(admin.ModelAdmin):
  list_display = ('name',)
  fields       = ('name', 'fivep', 'threep')
  search_fields = ('name', 'fivep', 'threep')

#############################################
@admin.register(Machine)
class MachineAdmin(admin.ModelAdmin):
  list_display = ('code', 'platform', 'name')
  fields       = ('code', 'platform', 'name')

#############################################
@admin.register(Peakcalling)
class PeakcallingAdmin(admin.ModelAdmin):
  # FIXME provenance here also.
  list_display    = ('code',
                     'factor_align_link', 'input_align_link')

  search_fields   = ('code',
                     'factor_align__lane__library__code',
                     'input_align__lane__library__code')
  readonly_fields = ('input_align', 'factor_align')

  def factor_align_link(self, obj):
    url = reverse('admin:osqpipe_alignment_change', args=(obj.factor_align.pk,))
    return '<a href="%s">%s</a>' % (url, obj.factor_align)
  factor_align_link.allow_tags = True
  
  def input_align_link(self, obj):
    url = reverse('admin:osqpipe_alignment_change', args=(obj.input_align.pk,))
    return '<a href="%s">%s</a>' % (url, obj.input_align)
  input_align_link.allow_tags = True

#############################################
@admin.register(Peakfile)
class PeakfileAdmin(admin.ModelAdmin):
  list_display    = ('filename', 'peakcalling_link', 'filetype')
  search_fields   = ('filename', 'filetype__name')
  readonly_fields = ('peakcalling',)
  fields = ('filename', 'peakcalling', 'filetype', 'checksum', 'description')

  def peakcalling_link(self, obj):
    url = reverse('admin:osqpipe_peakcalling_change', args=(obj.peakcalling.pk,))
    return '<a href="%s">%s</a>' % (url, obj.peakcalling)
  peakcalling_link.allow_tags = True

#############################################
@admin.register(Program)
class ProgramAdmin(admin.ModelAdmin):
  list_display  = ('program', 'version', 'type')
  fields        = ('program', 'version', 'type',
                   'options', 'files', 'description')

#############################################
class ProjectAdminForm(forms.ModelForm):

  '''A custom form allowing us to override the default height of
  multi-select boxes.'''

  class Meta:
    widgets = {
      'people': forms.SelectMultiple(attrs={ 'size': 15 })
      }

@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
  form          = ProjectAdminForm
  list_display  = ('name', 'code', 'description')
  fields        = ('name', 'code', 'lab', 'is_frozen',
                   'shortnames', 'filtered', 'description', 'people')

#############################################
@admin.register(QCfile)
class QCfileAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'laneqc_link', 'filetype')

  search_fields   = ('filename', 'filetype__name')
  readonly_fields = ('laneqc', 'date')
  
  fields = ('filename', 'checksum', 'filetype',
            'laneqc', 'date', 'description')

  def laneqc_link(self, obj):
    url = reverse('admin:osqpipe_laneqc_change', args=(obj.laneqc.pk,))
    return '<a href="%s">%s</a>' % (url, obj.laneqc)
  laneqc_link.allow_tags = True

#############################################
@admin.register(AlnQCfile)
class AlnQCfileAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'alignmentqc_link', 'filetype')

  search_fields   = ('filename', 'filetype__name')
  readonly_fields = ('alignmentqc', 'date')
  
  fields = ('filename', 'checksum', 'filetype',
            'alignmentqc', 'date', 'description')

  def alignmentqc_link(self, obj):
    url = reverse('admin:osqpipe_alignmentqc_change', args=(obj.alignmentqc.pk,))
    return '<a href="%s">%s</a>' % (url, obj.alignmentqc)
  alignmentqc_link.allow_tags = True

#############################################
@admin.register(HistologyImagefile)
class HistologyImagefileAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'sample_link', 'filetype')

  search_fields   = ('filename', 'filetype__name')
  readonly_fields = ('sample', 'date')
  
  fields = ('filename', 'checksum', 'filetype',
            'sample', 'date', 'block', 'batch', 'description')

  def sample_link(self, obj):
    url = reverse('admin:osqpipe_sample_change', args=(obj.sample.pk,))
    return '<a href="%s">%s</a>' % (url, obj.sample)
  sample_link.allow_tags = True

#############################################
@admin.register(Status)
class StatusAdmin(admin.ModelAdmin):
  list_display  = ('code', 'description', 'colour', 'authority')
  fields        = ('code', 'description', 'colour', 'authority',
                   'lanerelevant', 'libraryrelevant', 'sortcode')

#############################################
@admin.register(Strain)
class StrainAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name', 'description')

#############################################
@admin.register(Sex)
class SexAdmin(admin.ModelAdmin):
  list_display = ('__unicode__',)

#############################################
@admin.register(Tissue)
class TissueAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name', 'description')

#############################################
@admin.register(Condition)
class ConditionAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name', 'description')

#############################################
@admin.register(DoseUnit)
class DoseUnitAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name', 'description')

#############################################
@admin.register(SizeUnit)
class SizeUnitAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name', 'description')

#############################################
@admin.register(TreatmentAgent)
class TreatmentAgentAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description', 'accession')
  search_fields = ('name', 'description', 'accession')

#############################################
@admin.register(MergedAlignment)
class MergedAlignmentAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'genome')
  search_fields = ('alignment__lane__library__code',
                   'genome__code', 'alignment__lane__library__sample__name')
  readonly_fields = ('genome',)
  fields = ('genome',)

#############################################
@admin.register(MergedAlnfile)
class MergedAlnfileAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'alignment_link', 'filetype')

  search_fields   = ('filename', 'filetype__name')
  readonly_fields = ('alignment', 'date')
  
  fields = ('filename', 'checksum', 'filetype',
            'alignment', 'date', 'description')

  def alignment_link(self, obj):
    url = reverse('admin:osqpipe_mergedalignment_change', args=(obj.alignment.pk,))
    return '<a href="%s">%s</a>' % (url, obj.alignment)
  alignment_link.allow_tags = True

#############################################
@admin.register(ExternalRecord)
class ExternalRecordAdmin(admin.ModelAdmin):
  list_display = ('accession', 'repository', 'is_public')
  search_fields = ('accession', 'repository')

#############################################
@admin.register(ArchiveLocation)
class ArchiveLocationAdmin(admin.ModelAdmin):
  list_display  = ('name', 'host', 'root_path')
  search_fields = ('name', 'host', 'root_path')

  fields = ('name', 'host', 'root_path',
            'host_port', 'host_path', 'host_user',
            'host_delete_timelag')
