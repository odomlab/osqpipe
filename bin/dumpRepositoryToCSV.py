#!/usr/bin/env python

'''
Script to automatically dump a core set of library and lane annotation
out to CSV format for sharing with collaborators (via e.g. Dropbox).
'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqpipe.models import Lane, Library
from osqutil.config import Config

CONFIG = Config()

def lane_list_filetypes(lane):
  '''
  Return a string listing the available filetypes associated with the
  lane. This is somewhat stylised to only cover those file classes and
  types we think may be of interest.
  '''
  # I suppose we could implement these helpers within the model
  # classes themselves, but they're typically application-specific and
  # unlikely to be of general use.
  listing = []
  if any([ lf.filetype.code == 'fq'
           for lf in lane.lanefile_set.all() ]):
    listing += [u'fastq']
  if any([ alf.filetype.code == 'bam'
           for aln in lane.alignment_set.all()
           for alf in aln.alnfile_set.all() ]):
    listing += [u'bam']
  if any([ qcf.filetype.code == 'pdf'
           for lqc in lane.laneqc_set.all()
           for qcf in lqc.qcfile_set.all() ]):
    listing += [u'fastqc']

  return _helper_listing_to_string(listing)

def lane_aln_genomes(lane):
  '''
  Return a string listing the available genomes (and versions) used
  for alignments to this lane.
  '''
  listing = [ u"%s(%s)" % (aln.genome.code, aln.genome.version)
              for aln in lane.alignment_set.all() ]
  return _helper_listing_to_string(listing)

def lane_source_treatments(lane):
  '''
  Return a string listing the treatments applied to the source in this lane.
  '''
  listing = [ u"%s" % (treatment,)
              for treatment in lane.library.sample.source.sourcetreatment_set.all() ]
  return _helper_listing_to_string(listing)

def library_source_treatments(library):
  '''
  Return a string listing the treatments applied to the source in this library.
  '''
  listing = [ u"%s" % (treatment,)
              for treatment in library.sample.source.sourcetreatment_set.all() ]
  return _helper_listing_to_string(listing)

def lane_sample_characteristics(lane, category):
  '''
  Return a string listing the Characteristics with the specified
  category linked to the sample in this lane.
  '''
  listing = [ u"%s" % (char.value,)
              for char in lane.library.sample.characteristics.filter(category=category) ]
  return _helper_listing_to_string(listing)

def library_sample_characteristics(library, category):
  '''
  Return a string listing the Characteristics with the specified
  category linked to the sample in this library.
  '''
  listing = [ u"%s" % (char.value,)
              for char in library.sample.characteristics.filter(category=category) ]
  return _helper_listing_to_string(listing)

def library_projects(library, skipdefault=True):
  '''
  Return a string listing the projects to which this library
  belongs. Use skipdefault=False to include the default "all" project
  as specified in the main osqutil config file.
  '''
  prjs = library.projects.all()
  if skipdefault:
    prjs = prjs.exclude(code=CONFIG.defaultproject)
  listing = [ u"%s" % (project.name,) for project in prjs ]
  return _helper_listing_to_string(listing)

def lane_aln_programs(lane):
  '''
  Return a string listing the available programs (and versions) used
  for alignments to this lane.
  '''
  listing = []
  for aln in lane.alignment_set.all():
    listing += [ u"%s(%s)" % (prov.program.program, prov.program.version)
                 for prov in aln.provenance.all().order_by('rank_index') ]
  return _helper_listing_to_string(listing)

def lane_public_accessions(lane):
  '''
  Simply returns a string listing the public accessions linked to this
  lane.
  '''
  listing = [ rec.accession for rec in lane.public_records ]
  return _helper_listing_to_string(listing)

def lane_max_mapped_percent(lane, attr='mapped_percent'):
  '''
  Return the maximum percent mapped reads for a given lane.
  '''
  percs = [ getattr(aln, attr) for aln in lane.alignment_set.all() ]
  if len(percs) == 0:
    return u'NA'
  else:
    return unicode(max(percs))

def _helper_listing_to_string(listing):
  '''
  Quickly convert a list to a string, handling empty listings
  appropriately.
  '''
  if len(listing) == 0:
    listing = [u'NA']

  return u','.join(listing)

def _helper_optional_value(value, attr=None):
  if value is not None:
    if attr is not None:
      return unicode(getattr(value, attr))
    else:
      return unicode(value)
  else:
    return u'NA'

class RepositoryDumper(object):
  '''
  Base class defining some simple behaviours of RepositoryDumper objects.
  '''
  __slots__ = ('mapping', 'separator', 'objects')

  def __init__(self, separator="\t"):
    self.mapping   = []
    self.separator = separator
    self.objects   = []

  def header_string(self):
    '''
    Returns a string to be used as file header, ready for writing to
    output file.
    '''
    return self.separator.join([ elem[0] for elem in self.mapping ])

  def object_to_string(self, lane):
    '''
    Returns a string representing the passed lane object, ready for
    writing to output file.
    '''
    return self.separator.join([ elem[1](lane) for elem in self.mapping ])

  def format(self):
    '''
    Generator method which iterates through each dumped object and formats it for output.
    '''
    for obj in self.objects:
      yield self.object_to_string(obj) + u'\n'

class RepositoryLaneDumper(RepositoryDumper):
  '''
  Class which dumps a specific set of metadata from all libraries and
  lanes in the repository to an output file. The way in which we
  denormalise these various database tables is almost certainly
  specific to the intended application: allowing our collaborators to
  see roughly what data we have.
  '''

  def __init__(self, project=None, releaseworthy=False, *args, **kwargs):
    super(RepositoryLaneDumper, self).__init__(*args, **kwargs)
    
    # This mapping controls the output; the list contains 2-element
    # tuples which are (header_string, lambda function to retrieve row
    # contents).
    self.mapping = [
      ('Library',                lambda x: x.library.code),
      ('Experiment',             lambda x: _helper_optional_value(x.library.chipsample)),
      ('Facility',               lambda x: x.facility.code),
      ('Lane Number',            lambda x: str(x.lanenum)),

      ('Sex',                    lambda x: str(_helper_optional_value(x.library.sample.source.sex))),
      ('Tissue/Cell Line',       lambda x: x.library.sample.tissue.name),
      ('Strain',                 lambda x: _helper_optional_value(x.library.sample.source.strain, 'name')),
      ('Sample ID',              lambda x: _helper_optional_value(x.library.sample.name)),
      ('Individual',             lambda x: _helper_optional_value(x.library.sample.source.name)),
      ('Condition',              lambda x: _helper_optional_value(x.library.condition, 'name')),
      ('Treatments',             lambda x: lane_source_treatments(x)),
      ('Diagnosis',              lambda x: lane_sample_characteristics(x, 'Diagnosis')),
      ('Library Type',           lambda x: x.library.libtype.name),
      ('Library ChIP Factor',    lambda x: _helper_optional_value(x.library.factor, 'name')),
      ('Library ChIP Antibody',  lambda x: str(_helper_optional_value(x.library.antibody))),
      ('Sequencing Adapter 1',   lambda x: str(_helper_optional_value(x.library.adapter, 'sequence'))),
      ('Sequencing Adapter 2',   lambda x: str(_helper_optional_value(x.library.adapter2, 'sequence'))),

      ('Paired/Single Ended',    lambda x: 'PE' if x.paired else 'SE'),
      ('Flowcell ID',            lambda x: x.flowcell),
      ('Flowcell Lane',          lambda x: str(x.flowlane)),
      ('Run Date',               lambda x: str(x.rundate)),
      ('Read Length',            lambda x: str(x.readlength)),
      ('Total Reads Passed PF',  lambda x: str(x.total_passedpf)),

      # Consider spliting this into one filetype per column?
      ('Filetypes Available',    lambda x: lane_list_filetypes(x)),
      ('Alignment Genomes',      lambda x: lane_aln_genomes(x)),
      ('Alignment Programs',     lambda x: lane_aln_programs(x)),

      ('Max Percent Mapped',     lambda x: str(lane_max_mapped_percent(x))),
      ('Max Percent Uniquely Mapped', lambda x:\
         str(lane_max_mapped_percent(x, 'munique_percent'))),

      ('Marked As Failed',       lambda x: u'YES' if x.library.bad else u'no'),
      ('Comment',                lambda x: _helper_optional_value(x.library.comment)),
      ('Accessions',             lambda x: lane_public_accessions(x)),
    ]

    # This is broadly aping the approach in osqpipe.views to improve
    # performance and reduce overall load on the database. It's
    # maybe not perfect but it's a good place to start. Note that we
    # filter out virtual lanes (where passedpf=NULL).
    lanes = Lane.objects\
        .select_related('facility', 'library')\
        .prefetch_related('library__genome','library__libtype',
                          'library__factor','library__antibody',
                          'library__sample',
                          'library__sample__tissue',
                          'library__sample__source',
                          'library__sample__source__strain',
                          'library__sample__source__sex',
                          'library__adapter',
                          'lanefile_set', 'lanefile_set__filetype',
                          'laneqc_set', 'laneqc_set__qcfile_set',
                          'laneqc_set__qcfile_set__filetype',
                          'alignment_set', 'alignment_set__alnfile_set',
                          'alignment_set__alnfile_set__filetype',
                          'alignment_set__genome',
                          'alignment_set__provenance__program')\
        .exclude(passedpf__isnull=True)\
        .order_by('library__extra__code_text_prefix',
                  'library__extra__code_numeric_suffix')

    if project is not None:
      lanes = lanes.filter(library__projects__code=project)

    if releaseworthy is True:
      lanes = lanes.filter(library__release_worthy=True)

    self.objects = lanes

class RepositoryLibraryDumper(RepositoryDumper):
  '''
  Class which dumps a specific set of metadata from all libraries in
  the repository to an output file. The way in which we denormalise
  these various database tables is almost certainly specific to the
  intended application: allowing our collaborators to see roughly what
  data we have.

  '''

  def __init__(self, project=None, releaseworthy=False, *args, **kwargs):
    super(RepositoryLibraryDumper, self).__init__(*args, **kwargs)
    
    # This mapping controls the output; the list contains 2-element
    # tuples which are (header_string, lambda function to retrieve row
    # contents).
    self.mapping = [
      ('Library',                lambda x: x.code),
      ('Experiment',             lambda x: _helper_optional_value(x.chipsample)),

      ('Sex',                    lambda x: str(_helper_optional_value(x.sample.source.sex))),
      ('Tissue/Cell Line',       lambda x: x.sample.tissue.name),
      ('Strain',                 lambda x: _helper_optional_value(x.sample.source.strain, 'name')),
      ('Sample ID',              lambda x: _helper_optional_value(x.sample.name)),
      ('Individual',             lambda x: _helper_optional_value(x.sample.source.name)),
      ('Date of Birth',          lambda x: _helper_optional_value(x.sample.source.date_of_birth)),
      ('Date of Death',          lambda x: _helper_optional_value(x.sample.source.date_of_death)),
      ('Mother',                 lambda x: _helper_optional_value(x.sample.source.mother, 'name')),
      ('Father',                 lambda x: _helper_optional_value(x.sample.source.father, 'name')),
      ('Species',                lambda x: _helper_optional_value(x.sample.source.species, 'scientific_name')),
      ('Condition',              lambda x: _helper_optional_value(x.condition, 'name')),
      ('Treatments',             lambda x: library_source_treatments(x)),
      ('Diagnosis',              lambda x: library_sample_characteristics(x, 'Diagnosis')),
      ('TumourGrade',            lambda x: library_sample_characteristics(x, 'TumourGrade')),
      ('Library Type',           lambda x: x.libtype.name),
      ('Library ChIP Factor',    lambda x: _helper_optional_value(x.factor, 'name')),
      ('Library ChIP Antibody',  lambda x: str(_helper_optional_value(x.antibody))),
      ('Sequencing Adapter 1',   lambda x: str(_helper_optional_value(x.adapter, 'sequence'))),
      ('Sequencing Adapter 2',   lambda x: str(_helper_optional_value(x.adapter2, 'sequence'))),

      ('Marked As Failed',       lambda x: u'YES' if x.bad else u'no'),
      ('Projects',               lambda x: str(library_projects(x))),
      ('Comment',                lambda x: _helper_optional_value(x.comment)),
    ]

    # This is broadly aping the approach in osqpipe.views to improve
    # performance and reduce overall load on the database. It's
    # maybe not perfect but it's a good place to start. Note that we
    # filter out virtual lanes (where passedpf=NULL).
    libraries = Library.objects\
        .prefetch_related('genome','libtype',
                          'factor','antibody',
                          'sample',
                          'sample__tissue',
                          'sample__source',
                          'sample__source__strain',
                          'sample__source__sex',
                          'adapter',
                          'lane_set')\
        .exclude(lane__passedpf__isnull=True)\
        .distinct()\
        .order_by('extra__code_text_prefix',
                  'extra__code_numeric_suffix')

    if project is not None:
      libraries = libraries.filter(projects__code=project)

    if releaseworthy is True:
      libraries = libraries.filter(release_worthy=True)

    self.objects = libraries

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='Dump library metadata to CSV output.')

  PARSER.add_argument('-o', '--output', dest='output', type=str, required=True,
                      help='The name of the output file.')

  PARSER.add_argument('-p', '--project', dest='project', type=str, required=False,
                      help='A project code by which to filter the output.')

  PARSER.add_argument('--by-library', dest='bylib', action='store_true',
                      help='Dump data at the level of Library, rather than the default Lane level.')

  PARSER.add_argument('--release-only', dest='releaseworthy', action='store_true',
                      help='Only dump out data which has been linked to those libraries'
                      + 'specifically marked for release in the repository.')

  ARGS = PARSER.parse_args()

  if ARGS.bylib:
    DUMPER = RepositoryLibraryDumper(project=ARGS.project,
                                     releaseworthy=ARGS.releaseworthy)
  else:
    DUMPER = RepositoryLaneDumper(project=ARGS.project,
                                  releaseworthy=ARGS.releaseworthy)

  with open(ARGS.output, 'w') as outfh:
    outfh.write(DUMPER.header_string() + "\n")
    for lanestr in DUMPER.format():
      outfh.write(lanestr.encode('utf-8'))
    
