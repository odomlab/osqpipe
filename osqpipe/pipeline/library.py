#!/usr/bin/env python
#
# $Id$

'''Given a set of metadata on the command-line, create a new row in
the library repository table.'''

import sys
import re
import copy
from config import Config

from ..models import Factor, Genome, Antibody, Strain, Sex, Tissue, \
    Library, Libtype, Project, Adapter, Linkerset

from django.db import transaction

from setup_logs import configure_logging
LOGGER = configure_logging()

CONFIG = Config()

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

  def _add_if_confirmed(self, thing, thingcls):
    '''Add a controlled term to the database, if it's confirmed by the
    user.'''
    thingset = thingcls.objects.all()
    okay = self._confirm_field(thing,
                               thingset)
    if okay and not self.test_mode:
      thing.save()
    else:
      sys.exit("%s CV term not recognised as valid." %
               (thingcls.__name__))

  def _retrieve_cv(self, value, cls, user_key='controlled_name', **kwargs):
    '''Retrieve a database CV or create a new one if confirmed by the
    user. The default user key column is "controlled_name".'''
    if self.fuzzy:
      querykey = "%s__fuzzy" % user_key
    else:
      querykey = user_key
    kwargs[querykey] = value
    db_values = cls.objects.filter(**kwargs)
    if len(db_values) == 0:
      db_value = cls(**{ user_key: value })
      self._add_if_confirmed(db_value, cls)
      return(db_value)
    elif len(db_values) == 1:
      return db_values[0]
    else:
      raise ValueError("Multiple fuzzy CV matches found in database: %s"
                       % ( ", ".join([x.controlled_name for x in db_values])))

  @transaction.commit_on_success
  def _save_lib_to_database(self, code, keys, projects):
    '''
    Method to save library to the database, and link it with the
    appropriate projects. It is this project link which is the sole
    reason for wrapping this method in a database transaction.
    '''
    lib = Library(**keys)
    if not self.test_mode:
      LOGGER.info("Saving library to database: %s", code)
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
        raise StandardError("Requested project not found in the database: %s"
                            % proj)

    keys['code']   = code
    keys['genome'] = genome
    keys['tissue'] = tissue

    # Generate an ad-hoc description string.
    if 'description' not in keys or keys['description'] is None:
      fac = keys['factor'] if 'factor' in keys else 'unk'
      tis = keys['tissue'] if 'tissue' in keys else 'unk'
      ant = keys['antibody'] if 'antibody' in keys else 'unk'
      ind = keys['individual'] if 'individual' in keys else ''
      sta = keys['strain'] if 'strain' in keys else ''
      keys['description'] = ("%s_%s_%s_%s%s%s"
                             % (fac, tis, ant, keys['genome'], sta, ind))

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

      if 'adapter' in keys and keys['adapter']:
        keys['adapter'] = self._retrieve_cv(keys['adapter'], Adapter)

      if 'linkerset' in keys and keys['linkerset']:
        keys['linkerset'] = self._retrieve_cv(keys['linkerset'], Linkerset)

    except SystemExit, err:
      err = "Unable to load library %s: %s" % (code, err)
      if self.interactive:
        sys.exit(err)
      else:
        LOGGER.warning(err)
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

    # Write the library to the database in a managed transaction.
    lib = self._save_lib_to_database(code, keys, projects)

    return lib

