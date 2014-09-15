'''Configuration for the auto-generated django admin web interface.'''

from models import *
from django.contrib import admin
from django.core.urlresolvers import reverse
from django import forms

# Admin-level display config is held in each *Admin class; we register
# all model components except for LibraryNameIgnore which I suspect
# could be deprecated.

#############################################
class AdapterAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'sequence')
  search_fields = ('code', 'sequence')

admin.site.register(Adapter, AdapterAdmin)

#############################################
#FIXME we need to handle provenance (program/parameters) here
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
  
admin.site.register(Alignment, AlignmentAdmin)

#############################################
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

admin.site.register(Alnfile, AlnfileAdmin)

#############################################
class AntibodyAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'lot_number', 'description')
  search_fields = ('name', 'lot_number')

admin.site.register(Antibody, AntibodyAdmin)

#############################################
class FacilityAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'name', 'description')

admin.site.register(Facility, FacilityAdmin)

#############################################
class FactorAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name',)

admin.site.register(Factor, FactorAdmin)

#############################################
class FiletypeAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'name', 'suffix', 'description')

admin.site.register(Filetype, FiletypeAdmin)

#############################################
class GenomeAdmin(admin.ModelAdmin):
  list_display  = ('__unicode__', 'common_name', 'scientific_name', 'version')
  
  search_fields = ('code', 'common_name', 'scientific_name')

  fields = ('code', 'scientific_name', 'common_name', 'fasta',
            'fasta_md5sum', 'url', 'notes', 'version', 'blastdb')

admin.site.register(Genome, GenomeAdmin)

#############################################
class LaneAdmin(admin.ModelAdmin):
  list_display  = ('library', 'facility', 'lanenum',
                   'flowcell', 'runnumber', 'rundate')

  search_fields = ('library__code', 'facility__name', 'flowcell', 'runnumber')

  readonly_fields = ('external_records',)

  fieldsets = [
    (None,
     {'fields': ['library', 'facility', 'lanenum']}),
    ('Run Data',
     {'fields': ['flowcell', 'runnumber', 'machine', 'rundate', 'paired', 'readlength']}),
    ('LIMS Record',
     {'fields': ['genomicssampleid', 'usersampleid']}),
    ('QC Data',
     {'fields': ['reads', 'passedpf', 'mapped', 'failed', 'status', 'summaryurl']}),
    ('Other Info',
     {'fields': ['external_records', 'notes'] })
    ]

admin.site.register(Lane, LaneAdmin)

#############################################
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
  
admin.site.register(LaneQC, LaneQCAdmin)

#############################################
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
  
admin.site.register(Lanefile, LanefileAdmin)

#############################################
class LibfileAdmin(admin.ModelAdmin):
  list_display       = ('__unicode__', 'library_link', 'filetype')

  search_fields = ('filename', 'library__code', 'filetype__name')

  fields = ('filename', 'library', 'checksum',
            'filetype', 'description')

  def library_link(self, obj):
    url = reverse('admin:osqpipe_library_change', args=(obj.library.pk,))
    return '<a href="%s">%s</a>' % (url, obj.library)
  library_link.allow_tags = True
  
admin.site.register(Libfile, LibfileAdmin)

#############################################
class LibraryAdminForm(forms.ModelForm):

  '''A custom form allowing us to override the default height of
  multi-select boxes.'''

  class Meta:
    model = Library
    widgets = {
      'projects': forms.SelectMultiple(attrs={ 'size': 15 })
      }

class LibraryAdmin(admin.ModelAdmin):
  form = LibraryAdminForm
  
  list_display       = ('__unicode__', 'genome_link', 'libtype', 'chipsample', 'strain',
                        'individual', 'tissue', 'factor', 'antibody')

  search_fields  = ('code', 'genome__code', 'libtype__name', 'sex__name',
                    'strain__name', 'individual', 'chipsample',
                    'factor__name', 'tissue__name', 'antibody__name')

  fieldsets = [
    (None,
     {'fields': ['code'] }),
    ('Sample Annotation',
     {'fields': ['strain', 'individual', 'tissue', 'sex',
                 'factor', 'antibody'] }),
    ('Technical Info',
     {'fields': ['libtype', 'genome', 'paired',
                 'adapter', 'linkerset', 'barcode', 'bad']}),
    ('Project Info',
     {'fields': ['projects', 'chipsample', 'comment']})
    ]

  def genome_link(self, obj):
    url = reverse('admin:osqpipe_genome_change', args=(obj.genome.pk,))
    return '<a href="%s">%s</a>' % (url, obj.genome)
  genome_link.allow_tags = True
  
admin.site.register(Library, LibraryAdmin)

#############################################
class LibraryNameMapAdmin(admin.ModelAdmin):
  list_display  = ('limsname', 'libname')
  fields        = ('limsname', 'libname')
  search_fields = ('limsname', 'libname')

admin.site.register(LibraryNameMap, LibraryNameMapAdmin)

#############################################
class LibtypeAdmin(admin.ModelAdmin):
  list_display = ('name', 'code', 'description')
  fields       = ('name', 'code', 'description')

admin.site.register(Libtype, LibtypeAdmin)

#############################################
class LinkersetAdmin(admin.ModelAdmin):
  list_display = ('name',)
  fields       = ('name', 'fivep', 'threep')
  search_fields = ('name', 'fivep', 'threep')

admin.site.register(Linkerset, LinkersetAdmin)

#############################################
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

admin.site.register(Peakcalling, PeakcallingAdmin)

#############################################
class PeakfileAdmin(admin.ModelAdmin):
  list_display    = ('filename', 'peakcalling_link', 'filetype')
  search_fields   = ('filename', 'filetype__name')
  readonly_fields = ('peakcalling',)
  fields = ('filename', 'peakcalling', 'filetype', 'checksum', 'description')

  def peakcalling_link(self, obj):
    url = reverse('admin:osqpipe_peakcalling_change', args=(obj.peakcalling.pk,))
    return '<a href="%s">%s</a>' % (url, obj.peakcalling)
  peakcalling_link.allow_tags = True

admin.site.register(Peakfile, PeakfileAdmin)

#############################################
class ProgramAdmin(admin.ModelAdmin):
  list_display  = ('program', 'version', 'type')
  fields        = ('program', 'version', 'type',
                   'options', 'files', 'description')

admin.site.register(Program, ProgramAdmin)

#############################################
class ProjectAdminForm(forms.ModelForm):

  '''A custom form allowing us to override the default height of
  multi-select boxes.'''

  class Meta:
    model = Project
    widgets = {
      'people': forms.SelectMultiple(attrs={ 'size': 15 })
      }

class ProjectAdmin(admin.ModelAdmin):
  form          = ProjectAdminForm
  list_display  = ('name', 'code', 'description')
  fields        = ('name', 'code', 'lab',
                   'shortnames', 'filtered', 'description', 'people')

admin.site.register(Project, ProjectAdmin)

#############################################
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

admin.site.register(QCfile, QCfileAdmin)

#############################################
class StatusAdmin(admin.ModelAdmin):
  list_display  = ('code', 'description', 'colour', 'authority')
  fields        = ('code', 'description', 'colour', 'authority',
                   'lanerelevant', 'libraryrelevant', 'sortcode')

admin.site.register(Status, StatusAdmin)

#############################################
class StrainAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name', 'description')

admin.site.register(Strain, StrainAdmin)

#############################################
class SexAdmin(admin.ModelAdmin):
  list_display = ('__unicode__',)

admin.site.register(Sex, SexAdmin)

#############################################
class TissueAdmin(admin.ModelAdmin):
  list_display = ('__unicode__', 'description')
  search_fields = ('name', 'description')

admin.site.register(Tissue, TissueAdmin)

#############################################
class ExternalRecordAdmin(admin.ModelAdmin):
  list_display = ('accession', 'repository', 'is_public')
  search_fields = ('accession', 'repository')

admin.site.register(ExternalRecord, ExternalRecordAdmin)

