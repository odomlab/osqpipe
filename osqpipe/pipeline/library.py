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

'''Given a set of metadata on the command-line, create a new row in
the library repository table.'''

import sys
import re
import copy
from osqutil.config import Config

from ..models import Factor, Genome, Antibody, Strain, Sex, Tissue, \
    Library, Libtype, Project, Adapter, Linkerset, Sample, Source, \
    Condition

from django.db import transaction

from osqutil.setup_logs import configure_logging
LOGGER = configure_logging()

CONFIG = Config()

######################################################################
class AnnotationMismatchError(ValueError):
  '''
  Custom exception class used to signal a mismatch between input and database.
  '''
  pass

######################################################################
class LibraryHandler(object):

  '''Class encapsulating the creation of new library records.'''

  __slots__ = ('fuzzy', 'interactive', 'test_mode')

  def __init__(self, fuzzy=False, interactive=True, test_mode=False):
    self.fuzzy = fuzzy
    self.interactive = interactive
    self.test_mode   = test_mode

  def _confirm_value(self, name, value):
    '''Ask the user whether creating a new value in the repository is
    acceptable (calling context defines where this might be). Return
    True/False on success/failure.'''
    if not self.interactive:
      return False
    msg = "Confirm value of '%s' for field '%s' (y/n)?\n" % (value, name)
    sys.stderr.write(msg)
    response = sys.stdin.readline().strip().lower()
    okay = response == "y" or response == "yes"
    if okay:
      sys.stderr.write("Accepted.\n")
    else:
      sys.stderr.write("Not accepted.\n")
    return okay

  def _confirm_field(self, thing, collection):
    '''Ask the user to confirm the addition of a controlled term to
    the database. A list of pre-loaded values is displayed.'''
    if not self.interactive:
      return False
    msg = ("Confirm addition of '%s' to table '%s' (y/n)? "
           % (str(thing), type(thing)))
    sys.stderr.write(msg)
    response = sys.stdin.readline().strip().lower()
    okay = response == "y" or response == "yes"
    if okay:
      if len(collection) > 0:
        sys.stderr.write("Known %s values are: %s\n"
                         % (type(thing), " ".join([ x.controlled_name
                                                    for x in collection ])))
        sys.stderr.write(
          "Are you absolutely sure it's not one of these you mean?\n")
        response = sys.stdin.readline().strip().lower()
        okay = response == "y" or response == "yes"
    return okay

  def _add_if_confirmed(self, thing, thingcls, user_key='controlled_name'):
    '''Add a controlled term to the database, if it's confirmed by the
    user.'''
    thingset = thingcls.objects.all()
    okay = self._confirm_field(thing,
                               thingset)
    if okay and not self.test_mode:
      thing.save()
    else:
      sys.exit("%s CV term not recognised as valid (%s)." %
               (thingcls.__name__, getattr(thing, user_key)))

  def _retrieve_cv(self, value, cls, user_key='controlled_name', **kwargs):
    '''Retrieve a database CV or create a new one if confirmed by the
    user. The default user key column is "controlled_name".'''
    # Try to find an exact match first; fall back to fuzzy matching if
    # that's desired.
    kwargs[user_key] = value
    db_values = cls.objects.filter(**kwargs)

    if len(db_values) == 0 and self.fuzzy:
      querykey = "%s__fuzzy" % user_key
      del kwargs[user_key]
      kwargs[querykey] = value
      db_values = cls.objects.filter(**kwargs)

    if len(db_values) == 0:
      db_value = cls(**{ user_key: value })
      self._add_if_confirmed(db_value, cls, user_key)
      return(db_value)
    elif len(db_values) == 1:
      return db_values[0]
    else:
      raise ValueError("Multiple fuzzy CV matches found in database: %s"
                       % ( ", ".join([x.controlled_name for x in db_values])))

  @transaction.atomic
  def _save_lib_to_database(self, code, keys, projects):
    '''
    Method to save library to the database, and link it with the
    appropriate projects. It is this project link which is the sole
    reason for wrapping this method in a database transaction. Sample
    and Source objects will be created as necessary.
    '''
    sample_fields = ['tissue']
    source_fields = ['strain', 'sex']
    namefield     = 'individual'

    sourcekeys = dict( (k, v) for (k, v) in keys.iteritems() if k in source_fields )

    samplekeys = dict( (k, v) for (k, v) in keys.iteritems() if k in sample_fields )
    samplekeys['name'] = str(keys[namefield])

    try:

      # The Sample table has (name, tissue) as a unique key, which
      # greatly simplifies the business of pulling out a correctly
      # annotated Sample object. We could remove sample.name
      # altogether and key by (source, tissue); however this would not
      # allow the HCC project special case to work as well.
      sample = Sample.objects.get(**samplekeys)

      # A quick check on the linked Source is advisable nonetheless.
      for field in source_fields:
        if field in sourcekeys:
          if getattr(sample.source, field) != sourcekeys[field]:
            raise AnnotationMismatchError("Probable mislabeled sample ID, Source fields"
                                          + " disagree with database: %s (%s)" % (sample.source.name, field))
        elif getattr(sample.source, field) is not None:

          # We have no annotation in sheet but annotation in
          # DB. Possibly we could ignore this error and just assume
          # that data has been omitted for brevity?
          raise AnnotationMismatchError("Probable mislabeled sample ID, Source fields"
                                        + " disagree with database: %s (%s)" % (sample.source.name, field))
      
    except Sample.DoesNotExist:

      # This now needs to also handle Source retrieval/creation:
      try:
        source = Source.objects.get(name=str(keys[namefield]))
        for field in source_fields:
          if getattr(source, field) != sourcekeys[field]:
            raise AnnotationMismatchError("Probable mislabeled sample ID, Source fields"
                                          + " disagree with database: %s (%s)" % (source.name, field))

      except Source.DoesNotExist:
        sourcekeys['name'] = str(keys[namefield])
        source = Source(**sourcekeys)
        if not self.test_mode:
          LOGGER.info("Saving source to database: %s", source.name)
          source.save()

      sample = Sample(**samplekeys)
      sample.source = source
    
    keys   = dict( (k, v) for (k, v) in keys.iteritems()
                   if k not in sample_fields + source_fields + [namefield] )
    lib    = Library(**keys)
    
    if not self.test_mode:
      LOGGER.info("Saving sample to database: %s", sample.name)
      sample.save()
      LOGGER.info("Saving library to database: %s", code)
      lib.sample = sample
      lib.save()

      # Add to the specified project(s).
      for proj in projects:
        lib.projects.add(proj)

    else:
      LOGGER.info("TEST MODE: Would have saved library to database: %s", code)

  def add(self, libtype, code, genome, tissue, projcodes=None, opts=None):

    '''Method takes a set of required metadata and a dict of optional
    metadata (opts).'''

    if projcodes is None:
      projcodes = []

    if opts is None:
      opts = {}

    keys = copy.copy(opts) # No unpleasant side effects in the calling code.

    # We add to the specified project(s), but also to an over-arching
    # project which will represent all libraries.
    defproj = CONFIG.defaultproject
    if defproj not in projcodes:
      projcodes.append( defproj )

    # Case-insenstive 'in' query is a bit clunky in Django 1.5 (see
    # http://stackoverflow.com/a/9464142):
    projects  = Project.objects.filter(
      code__iregex=r'(' + '|'.join([re.escape(p) for p in projcodes]) + ')')
    retrieved = [ x.code.lower() for x in projects ]
    if defproj.lower() not in retrieved:
      raise StandardError("Default project not found in the database: %s"
                          % defproj)
    for proj in projcodes:
      if proj.lower() not in retrieved:
        if self.interactive:
          raise StandardError("Requested project not found in the database: %s"
                              % proj)
        else:
          LOGGER.warning("Requested project not found in the database: '%s'. Skipping", proj)
          return

    keys['code']   = code
    keys['genome'] = genome
    keys['tissue'] = tissue

    print "Adapter=\"%s\"" % keys['adapter']
    
    # A sensible fallback to make sure samples are treated
    # appropriately.
    indivkey = 'individual'
    if indivkey not in keys or keys[indivkey] is None or keys[indivkey] == '':
      keys[indivkey] = code

    # If we're using approximate matching, check our synonyms
    # list. Note that we *could* store these directly in the
    # repository. FIXME?
    if self.fuzzy:
      genome = CONFIG.genome_synonyms.get(genome, genome)

    # Confirm an exact match for the genome (case insensitive).
    try:
      db_genome = Genome.objects.get(code__iexact=genome)
    except Genome.DoesNotExist, err:
      if self.interactive:
        sys.exit(
          ("ERROR: genome '%s' not in database (hint:"
           + " use the code, e.g. mmu, hsa-hg19).")
          % genome)
      else:
        LOGGER.warning("Genome for %s not found in database: '%s'. Skipping.",
                       code, genome)
        return

    keys['genome'] = db_genome
    
    # We allow configurably fuzzy matching for tissue, factor, antibody, strain.
    try:
      keys['tissue'] = self._retrieve_cv(keys['tissue'], Tissue)
      keys['libtype'] = self._retrieve_cv(libtype, Libtype, 'code')

      if 'factor' in keys and keys['factor']:
        keys['factor'] = self._retrieve_cv(keys['factor'], Factor)

      if 'antibody' in keys and keys['antibody']:
        if 'lot_number' not in keys:
          keys['lot_number'] = 'unknown'

        # This will often write to the database outside of a
        # transaction; I don't think that's a problem, though (if an
        # antibody lot exists, it exists, and may as well be
        # represented in the db; erroneous lot numbers are outside the
        # scope of this class).
        keys['antibody'] = self._retrieve_cv(keys['antibody'], Antibody,
                                            lot_number=keys['lot_number'],
                                            autocreate_lots=True,
                                            test_mode=self.test_mode)

      if 'strain' in keys and keys['strain']:
        keys['strain'] = self._retrieve_cv(keys['strain'], Strain)

      if 'sex' in keys and keys['sex']:
        keys['sex'] = self._retrieve_cv(keys['sex'], Sex)

      if 'condition' in keys and keys['condition']:
        keys['condition'] = self._retrieve_cv(keys['condition'], Condition)

      if 'adapter' in keys and keys['adapter']:
        keys['adapter'] = self._retrieve_cv(keys['adapter'], Adapter)

      if 'linkerset' in keys and keys['linkerset']:
        keys['linkerset'] = self._retrieve_cv(keys['linkerset'], Linkerset)

      if 'adapter2' in keys and keys['adapter2']:
        keys['adapter2'] = self._retrieve_cv(keys['adapter2'], Adapter)

    except SystemExit, err:
      err = "Unable to load library %s: %s" % (code, err)
      if self.interactive:
        sys.exit(err)
      else:
        LOGGER.error(err)
        return

    # Special case for antibody lot number.
    if 'lot_number' in keys:
      del(keys['lot_number'])

    # If library already present, skip it. Otherwise ask the user if
    # the code doesn't conform to our standard pattern.
    if Library.objects.filter(code__iexact=code).count() > 0:
      sys.exit("Code '%s' is already present in the database." % code)
    code_pat = re.compile(r"do\d{2,}")
    matchobj = code_pat.match(code)
    if not matchobj:
      okay = self._confirm_value("name", code)
      if not okay:
        if self.interactive:
          sys.exit("User terminated script.")
        else:
          LOGGER.warning("Unable to save library: unexpected code format: %s",
                         code)
          return

    # Write the sample and library to the database in a managed
    # transaction. This method knows which field goes where.
    lib = self._save_lib_to_database(code, keys, projects)

    return lib

