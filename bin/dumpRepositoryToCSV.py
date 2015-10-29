#!/usr/bin/env python

'''
Script to automatically dump a core set of library and lane annotation
out to CSV format for sharing with collaborators (via e.g. Dropbox).
'''

from osqpipe.models import Lane

from osqpipe.pipeline.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

def helper_list_filetypes(lane):
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
    listing += ['fastq']
  if any([ alf.filetype.code == 'bam'
           for aln in lane.alignment_set.all()
           for alf in aln.alnfile_set.all() ]):
    listing += ['bam']
  if any([ qcf.filetype.code == 'pdf'
           for lqc in lane.laneqc_set.all()
           for qcf in lqc.qcfile_set.all() ]):
    listing += ['fastqc']

  return _helper_listing_to_string(listing)

def helper_aln_genomes(lane):
  '''
  Return a string listing the available genomes (and versions) used
  for alignments to this lane.
  '''
  listing = [ "%s(%s)" % (aln.genome.code, aln.genome.version)
              for aln in lane.alignment_set.all() ]
  return _helper_listing_to_string(listing)

def helper_aln_programs(lane):
  '''
  Return a string listing the available programs (and versions) used
  for alignments to this lane.
  '''
  listing = []
  for aln in lane.alignment_set.all():
    listing += [ "%s(%s)" % (prov.program.program, prov.program.version)
                 for prov in aln.provenance.all().order_by('rank_index') ]
  return _helper_listing_to_string(listing)

def helper_public_accessions(lane):
  '''
  Simply returns a string listing the public accessions linked to this
  lane.
  '''
  listing = [ rec.accession for rec in lane.public_records ]
  return _helper_listing_to_string(listing)

def helper_max_mapped_percent(lane, attr='mapped_percent'):
  '''
  Return the maximum percent mapped reads for a given lane.
  '''
  percs = [ getattr(aln, attr) for aln in lane.alignment_set.all() ]
  if len(percs) == 0:
    return 'NA'
  else:
    return max(percs)

def _helper_listing_to_string(listing):
  '''
  Quickly convert a list to a string, handling empty listings
  appropriately.
  '''
  if len(listing) == 0:
    listing = ['NA']

  return ",".join(listing)

def _helper_optional_value(value, attr=None):
  if value is not None:
    if attr is not None:
      return getattr(value, attr)
    else:
      return value
  else:
    return 'NA'

class RepositoryDumper(object):
  '''
  Class which dumps a specific set of metadata from all libraries and
  lanes in the repository to an output file. The way in which we
  denormalise these various database tables is almost certainly
  specific to the intended application: allowing our collaborators to
  see roughly what data we have.
  '''
  __slots__ = ('mapping', 'separator')

  def __init__(self, separator="\t"):
    # This mapping controls the output; the list contains 2-element
    # tuples which are (header_string, lambda function to retrieve row
    # contents).
    self.mapping = [
      ('Library',                lambda x: x.library.code),
      ('Facility',               lambda x: x.facility.code),
      ('Lane Number',            lambda x: str(x.lanenum)),

      ('Sex',                    lambda x: str(_helper_optional_value(x.library.sex))),
      ('Tissue/Cell Line',       lambda x: x.library.tissue.name),
      ('Strain',                 lambda x: _helper_optional_value(x.library.strain, 'name')),
      ('Library Type',           lambda x: x.library.libtype.name),
      ('Library ChIP Factor',    lambda x: _helper_optional_value(x.library.factor, 'name')),
      ('Library ChIP Antibody',  lambda x: str(_helper_optional_value(x.library.antibody))),

      ('Paired/Single Ended',    lambda x: 'PE' if x.paired else 'SE'),
      ('Flowcell ID',            lambda x: x.flowcell),
      ('Flowcell Lane',          lambda x: str(x.flowlane)),
      ('Run Date',               lambda x: str(x.rundate)),
      ('Read Length',            lambda x: str(x.readlength)),
      ('Total Reads Passed PF',  lambda x: str(x.total_passedpf)),

      # Consider spliting this into one filetype per column?
      ('Filetypes Available',    lambda x: helper_list_filetypes(x)),
      ('Alignment Genomes',      lambda x: helper_aln_genomes(x)),
      ('Alignment Programs',     lambda x: helper_aln_programs(x)),

      ('Max Percent Mapped',     lambda x: str(helper_max_mapped_percent(x))),
      ('Max Percent Uniquely Mapped', lambda x:\
         str(helper_max_mapped_percent(x, 'munique_percent'))),

      ('Marked As Failed',       lambda x: 'YES' if x.library.bad else 'no'),
      ('Comment',                lambda x: _helper_optional_value(x.library.comment)),
      ('Accessions',             lambda x: helper_public_accessions(x)),
      ]
    self.separator = separator

  def header_string(self):
    '''
    Returns a string to be used as file header, ready for writing to
    output file.
    '''
    return self.separator.join([ elem[0] for elem in self.mapping ])

  def laneobj_to_string(self, lane):
    '''
    Returns a string representing the passed lane object, ready for
    writing to output file.
    '''
    return self.separator.join([ elem[1](lane) for elem in self.mapping ])

  def dump_to_file(self, outfile):
    '''
    Dumps everything to the specified output file.
    '''
    with open(outfile, 'w') as outfh:
      outfh.write(self.header_string() + "\n")

      # This is broadly aping the approach in osqpipe.views to improve
      # performance and reduce overall load on the database. It's
      # maybe not perfect but it's a good place to start. Note that we
      # filter out virtual lanes (where passedpf=NULL).
      lanes = Lane.objects\
          .select_related('facility', 'library')\
          .prefetch_related('library__genome', 'library__tissue','library__libtype',
                            'library__strain','library__factor','library__antibody',
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
      for lane in lanes:
        lanestr = self.laneobj_to_string(lane)
        outfh.write(lanestr + "\n")

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='Dump library metadata to CSV output.')

  PARSER.add_argument('-o', '--output', dest='output', type=str, required=True,
                      help='The name of the output file.')

  ARGS = PARSER.parse_args()

  DUMPER = RepositoryDumper()

  DUMPER.dump_to_file(ARGS.output)
