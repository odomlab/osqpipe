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

'''Custom model managers.'''

####################################################################################
# ControlledVocab management. Primarily concerned with handling
# phonetic ("fuzzy") searches. We only support fuzzy searching on
# filter, not get (since a fuzzy search may return multiple hits
# anyway).

from django.db import models
from django.db.models import Max
import re
import os

from osqutil.utilities import get_filename_libcode

from osqutil.setup_logs import configure_logging
LOGGER = configure_logging('models')

from osqutil.config import Config
CONFIG = Config()

class ControlledVocabManager(models.Manager):

  # This will become a dict of dicts. N.B. deactivated as the caching mechanism currently creates more problems than it solves.
#  _fuzzyCache = {}
  _fuzzy_re    = re.compile('(.*)__fuzzy$')

  def filter(self, **kwargs):

    kwargs = self._map_controlled_field(kwargs)
    (fuzzy_args, filter_args) = self.extract_args_regex(kwargs, self._fuzzy_re)
    objs = super(ControlledVocabManager, self).filter(**filter_args)

    if len(fuzzy_args) > 0:
      objs = self._fuzzy_filter(fuzzy_args, objs)

    return objs

  def get(self, **kwargs):

    kwargs = self._map_controlled_field(kwargs)
    return super(ControlledVocabManager, self).get(**kwargs)

  def _map_controlled_field(self, kwargs):
    '''Replace 'controlled_name' in arguments with the actual model
    field to be queried.'''
    cleaned = {}
    cont_re = re.compile('^controlled_name(.*)')
    for (k, v) in kwargs.iteritems():
      cont_match = cont_re.match(k)
      if cont_match:
        field = self.model._controlled_field
        if field is None:
          raise StandardError("Controlled field not designated for class %s"
                              % self.model)
        else:
          cleaned[ field + cont_match.group(1) ] = v
      else:
        cleaned[k] = v
    return cleaned

  def extract_args_regex(self, kwargs, regex):
    '''Split our argument list into those matching a regex, and
    everything else.'''
    matching_args  = {}
    other_args     = {}

    for (k, v) in kwargs.iteritems():
      matching = regex.match(k)
      if matching:
        matching_args[ matching.group(1) ] = v
      else:
        other_args[ k ] = v
    return(matching_args, other_args)

  def _fuzzy_filter(self, fuzzy_args, objs):
    '''Query the database using a phonetic matching scheme. Currently
    marked as internal-use-only as this is potentially quite expensive.'''

    if objs.count() == 0 or len(fuzzy_args) == 0:
      return objs

    # We could also use soundex or nysiis here. I've chosen Double
    # Metaphone because it's worked well in similar projects.
    import fuzzy
    dmeta = fuzzy.DMetaphone()

    cname = objs[0].__class__.__name__
    cache = {} # self.__class__._fuzzyCache

    # Set up the class-specific sub-dict.
    if cname not in cache:
      cache[cname] = {}

    # Iterate over all returned objects, cache phonetic terms for
    # each. This is the major source of inefficiency in this method.
    for obj in objs:
      for attr in fuzzy_args.iterkeys():
        attrhash = cache[cname].setdefault(attr, {})
        if obj not in attrhash:
          fuzz = self._fuzzy_encode(vars(obj)[attr], dmeta)
          attrhash[obj] = fuzz

    # Reverse the mapping hash.
    revmap = {} # A dict (by attr) of dicts (by fuzz) of sets (of dbids)
    for (attr, attrhash) in cache[cname].iteritems():
      revhash = revmap.setdefault(attr, {})
      for (dbid, fuzz) in attrhash.iteritems():
        objfuzz = revhash.setdefault(fuzz, set())
        objfuzz.add(dbid)

    # Match our user key to the phonetic terms from the database.
    members = []
    for (attr, value) in fuzzy_args.iteritems():
      fuzz = self._fuzzy_encode(value, dmeta)
      if attr in revmap:
        if fuzz in revmap[attr]:
          hits = revmap[attr][fuzz]
          LOGGER.info("Matched requested %s %s=%s to %s in the database.",
                      cname, attr, value,
                      ";".join([ vars(x)[attr] for x in hits ]))
          members += list(hits)

    # Filter out any objects deleted since loading the cache. This
    # also has the effect of making members a proper queryset again.
    if len(members) > 0:
      members = self.get_queryset().filter( pk__in=[ x.pk for x in members ] )
    else:
      members = self.none()

    # Detect potential collisions here, at the point of output.
    uniq    = set( x.controlled_name for x in members )
    membstr = ", ".join(uniq)
    if len(uniq) > 1:
      LOGGER.warning("Possible collision in %s phonetic representations:"
                     + " %s all represented as %s.",
                     cname, membstr, fuzz)

    return members
        
  @staticmethod
  def _fuzzy_encode(string, dmeta, minlen=6):
    '''Extension to whatever phonetic scheme we're using to explicitly
    capture digits; this seems to be particularly important for our
    use-case and without it we get far too many collisions.'''

    string = unicode(string) # because Excel can't constrain its variable types

    # Short strings really aren't helped by this approach. We just
    # strip out the most common indel characters.
    if len(string) <= minlen:
      return re.sub(r'[-_ /\\]', '', string).lower()

    fuzz = ''
    parts = re.findall(r'[A-Za-z]+|\d+', string)
    for token in parts:
      if re.match(r'\d', token):
        fuzz = fuzz + token
      elif len(token) > 0:
        phon = dmeta(token)[0]
        fuzz = fuzz + (phon if phon is not None else token)
    return fuzz.lower() # case insensitive

class AntibodyManager(ControlledVocabManager):

  _lotnum_re = re.compile('^(lot_number(?!__fuzzy).*)$')

  def get(self, autocreate_lots=False, test_mode=False, **kwargs):

    (lotnum_args, _stripped) = self.extract_args_regex(kwargs, self._lotnum_re)
    if len(lotnum_args) == 0:
      kwargs['lot_number'] = 'unknown'
      return super(AntibodyManager, self).get(**kwargs)
    else:
      objs = self.filter(autocreate_lots=autocreate_lots,
                         test_mode=test_mode,
                         **kwargs)
      num = objs.count()
      if num == 0:
        raise self.model.DoesNotExist("No antibody found for query: %s" % (kwargs))
      elif num > 1:
        raise self.model.MultipleObjectsReturned(
          "get() returned more than one %s -- it returned %s! Lookup parameters were %s"
          % (self.model._meta.object_name, num, kwargs))
      else:
        return objs[0]

  def filter(self, autocreate_lots=False, test_mode=False, **kwargs):

    (lotnum_args, kwargs) = self.extract_args_regex(kwargs, self._lotnum_re)
    objs = super(AntibodyManager, self).filter(**kwargs)

    # If we're fuzzy-matching on lot number it's unlikely we're going
    # to want to support auto-creation (I hope). I guess we may
    # revisit this FIXME.
    if 'lot_number__fuzzy' in kwargs:
      if autocreate_lots:
        LOGGER.warn("Autocreation of lot numbers not supported"
                    + " when using fuzzy queries on lot numbers.")
      return objs

    if len(lotnum_args) == 0:
      lotnum_args['lot_number'] = 'unknown' # default lot_number

    antibodies = objs.filter(**lotnum_args)

    if antibodies.count() == 0:

      # "standard" queries failed.
      if objs.count() > 0 and autocreate_lots and 'lot_number' in lotnum_args:

        # This is more of a safety measure; the alternative would be too powerful.
        names = set([ x.name for x in objs ])
        if len(names) > 1:
          LOGGER.error("Autocreation of lot numbers not supported"
                       + " on multiple antibodies at once. Be more specific.")
          return objs.none()

        lotnum  = lotnum_args['lot_number']
        new_obj = self.model(name=objs[0].name, lot_number=lotnum)
        message = ("new lot_number for Antibody %s in database: %s"
                   % (objs[0].name, lotnum))

        if test_mode:
          LOGGER.warn("Would have created %s", message)
        else:
          LOGGER.warn("Creating %s", message)
          new_obj.save()

          # Re-query seems to be the best way to retrieve the new
          # object in a queryset.
          antibodies = super(AntibodyManager, self).filter(lot_number=lotnum, **kwargs)

    return antibodies

####################################################################################
# Custom Filetype manager class.

class FiletypeManager(models.Manager):

  def guess_type(self, fname):

    '''Attempt te match the suffix of the supplied filename with the
    file types represented in the database.'''

    fnparts = os.path.splitext(fname)

    if fnparts[1] == CONFIG.gzsuffix:
      realfn = fnparts[0]
    else:
      realfn = fname

    fn_suff = os.path.splitext(realfn)[1]

    return self.get(suffix=fn_suff)

####################################################################################
# Custom Library manager class.

class LibraryManager(models.Manager):

  def search_by_name(self, name):

    '''Fuzzy matching approach to try and retrieve a library code from
    a pseudo-structured string.'''

    from models import LibraryNameMap

    # basic policy: try to match whole name (possibly removing suffix and
    # lane ID ("CRI01" etc)).  Then try to match just the stuff preceding the
    # first underscore.  Then try looking up both in the library_name_map.

    base = os.path.splitext(name)[0]
    LOGGER.debug("search_by_name: base='%s'", base)

    lanePat = re.compile(r"(.*)_(\w{3}\d{2})$")
    matchobj = lanePat.match(base)
    if matchobj:
      libname = matchobj.group(1)
      lane = matchobj.group(2)
      LOGGER.debug("found lane code: '%s' -> '%s'", base, lane)
    else:
      libname = name

    lib = None

    try:
      lib = self.get(code__iexact=libname)
      LOGGER.debug("search_by_name: '%s' -> '%s'", name, lib.code)

    except self.model.DoesNotExist, err:
      pref = libname.split("_")[0]

      try:
        lib = self.get(code__iexact=pref)
        LOGGER.debug("search_by_name: '%s' -> '%s' -> '%s'",
                      name, pref, lib.code)

      except self.model.DoesNotExist, err:

        # FIXME it would be nice to deprecate this at some point.
        try:
          mapping = LibraryNameMap.objects.get(limsname=libname)
          lib = self.get(code__iexact=mapping.libname)
          if lib:
            LOGGER.debug("search_by_name: '%s' -> '%s' -> '%s'",
                          name, libname, lib.code)

        except LibraryNameMap.DoesNotExist, err:
          pass

    if lib is None:
      raise self.model.DoesNotExist("Library not found for name: %s" % name)
    
    return lib

  def search_by_filename(self, fname):
    '''Retrieve a Library object for a given filename from the
    repository database. This will cycle through the current expected
    naming scheme, an alternate (older?) scheme, and finally tries to
    match based on the Solexa LIMS naming scheme.'''

    from models import LibraryNameMap

    pref = get_filename_libcode(fname)  # Current naming scheme.
    suff = os.path.splitext(fname)[1]
    lib  = None

    try:
      lib  = self.get(code=pref)

    except self.model.DoesNotExist, err: # Unsure? old scheme?
      pat = re.compile(r"(.*)_(\w{3}\d{2})%s" % (suff, ))
      matchobj = pat.match(os.path.basename(fname))

      try:
        if matchobj:
          base = matchobj.group(1)
          lib = self.get(code=base)

      except self.model.DoesNotExist, err:

        # Solexa LIMS
        try:
          mapping = LibraryNameMap.objects.get(limsname=pref)
          lib     = self.get(code=mapping.libname)

        except LibraryNameMap.DoesNotExist, err:
          pass # But throw exception for a missing Library!

    if lib is None:
      raise self.model.DoesNotExist("Library not found for filename: %s" % fname)
    
    return lib

# Keeping this global allows for all the lane numbers assigned within
# a given process to avoid collisions. To avoid collisions with other
# processes, use transactions.
LIBRARY_LANE_CACHE = {}

class LaneManager(models.Manager):

  def next_lane_number(self, library):
    '''Return the next available flowcell lane number for a given
    library, starting at 1. Please note that new lane numbers are not
    stored in the database until any newly-created Lanes are saved,
    and therefore the only safe way to use this method is from within
    a transaction which retrieves the new number and saves the new
    Lane object.'''
    if library in LIBRARY_LANE_CACHE:
      currTop = LIBRARY_LANE_CACHE[library]
    else:
      currTop = super(LaneManager, self).filter(library=library).aggregate(Max('lanenum'))['lanenum__max']
    if currTop is None:
      currTop = 0
    newTop = currTop + 1
    LIBRARY_LANE_CACHE[library] = newTop
    return newTop

