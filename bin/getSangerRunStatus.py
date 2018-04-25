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

# Basic design:
#
# First, we set up an ssh tunnel to the Sanger squid proxy, using
# sshSangerTunnel.py
#
# We query http://npg.sanger.ac.uk/perl/npg/runs
# using our project name and pull down run lanes, status, date for all
# status types; munge the resulting table (we need to know the defined
# order of Sanger status to figure out the newest for a given
# lane). Store the lane alongside its latest status in our database.
#
# Failure mode: if the web site can't be accessed, email admins as
# it's likely the ssh tunnel has gone down.
#
# Design feature: every status has an "authority" field
# (i.e. 'Sanger', 'CRI Genomics', or whatever). This is just a link
# into the Facility table. Each authority has a "trigger" status at
# which point it hands off to the internal, 'local' (null) authority
# (e.g.: 'Sanger:run archived' or possibly 'Sanger:qc complete'). Each
# authority is effectively handed control over its own status flags
# (we may need to cook some up for CRI; see what the LIMS has).
# 
# Note that a run can breeze through the trigger status and get to the
# next one, so this system includes a DB field ('sortcode') containing
# status rank which should help to detect status downstream of the
# triggers.

import os
import re
import requests
from urlparse import urljoin, urlparse
from bs4 import BeautifulSoup
from datetime import datetime
from lxml import etree as ET

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from django.db import transaction
from osqpipe.models import Lane, Facility, Library, Status, Machine
from osqutil.config import Config
from osqpipe.pipeline.smtp import email_admins

# Script-specific globals; no real need to put these in the main config.
CACHE_FILE = os.path.join(os.path.expanduser('~'), '.sanger_run_status_latest')

# This needs to match the format stored in the CACHE_FILE and on the
# Sanger web site itself.
DATEFORMAT = '%Y-%m-%d %H:%M:%S'

# Number of seconds to wait on the web requests.
REQ_TIMEOUT = 300

# This is the link queried by embedded AJAX in the target web pages;
# see the output of Firebug (Firefox extension) to see how it's used
# (and if it's been changed).
#SANGER_INTAPI_URL = 'http://psd-support.internal.sanger.ac.uk:6600/'

class SangerLimsManager(object):
  '''Context manager class for the SangerLims worker class.'''

  def __init__(self, *args, **kwargs):
    self.lims = SangerLims(*args, **kwargs)

  def __enter__(self):
    return self.lims

  def __exit__(self, exc_type, exc_value, traceback):
    '''Write self.lims.running_status_date to cache file, if no
    exception has been raised.'''
    if exc_type is None:
      LOGGER.info("Caching latest status date.")
      with open(self.lims.cache_file, 'w') as cache:
        cache.write(self.lims.running_status_date.strftime(DATEFORMAT))

class Table:
  '''Simple iterator class.'''
  def __init__(self, table):
    self.rows = [x.split("\t") for x in table.splitlines()]
    self.current = 0

  def __iter__(self):
    return self

  def next(self):
    if self.current >= len(self.rows):
      raise StopIteration
    else:
      self.current += 1
      return self.rows[self.current - 1]
    
class HtmlTable(Table):
  '''Simple iterator class to make handling HTML tables slightly simpler.'''
  def __init__(self, table):
    self.rows = table.find_all('tr')
    self.current = 0

class SangerLims(object):
  '''The primary worker class used to run queries of the Sanger LIMS
  website.'''
  __slots__ = ('last_status_date', 'running_status_date', 'cache_file', 'conf', '_study_cache', '_missing_libcodes')

  def __init__(self):

    self.cache_file = CACHE_FILE
    self.conf = Config()
    self._study_cache = {}
    self._missing_libcodes = set()

    # This defines the date from which all rows will be checked for
    # new information.
    if os.path.exists(self.cache_file):
      with open(self.cache_file, 'r') as cache:
        date = cache.readline().strip()
        self.last_status_date = datetime.strptime(date, DATEFORMAT)
    else:
      # Fallback if cache file not present.
      self.last_status_date = datetime.fromtimestamp(0)

    # The running_status_date just keeps track of the most recent
    # status_date seen. It will be put into the cache file upon exit.
    self.running_status_date = datetime.fromtimestamp(0)

  def _parse_param_href(self, href, param):
    '''Parse an anchor element to pull out the ID for a given
    parameter and return it in a tuple alongside the displayed anchor
    text.'''
    match = re.search(r'%s=(\d+|all)' % param, href['href'])
    if match:
      return((href.text, match.group(1)))
    else:
      raise StandardError("Unable to parse run status href URL.")

  def _run_status_mapping(self, soup):
    '''Build a run status to param ID mapping dict.'''
    runspan = soup.find('span', text='Runs')
    hrefs   = runspan.parent.find_all('a')
    mapping = dict( self._parse_param_href(a, 'id_run_status_dict') for a in hrefs )
    mapping = dict( (k, int(v)) for (k, v) in mapping.iteritems() if k != 'all' )
    return(mapping)

  def _instrument_mapping(self, soup):
    '''Build an instrument type to param ID mapping dict.'''
    instspan = soup.find('span', text='Instrument Types')
    hrefs    = instspan.parent.find_all('a')
    mapping  = dict( self._parse_param_href(a, 'id_instrument_format') for a in hrefs )
    mapping  = dict( (k, int(v)) for (k, v) in mapping.iteritems() if k != 'all' )
    return(mapping)

  # see http://stackoverflow.com/a/3668771
  def _get_content(self, url, *args, **kwargs):
    '''Wrapper method for requests.get which adds in http proxy
    support and automatically follows HTML redirect tags.'''
    proxies = { 'http': self.conf.sanger_http_proxy }
    LOGGER.debug("Querying url %s", url)
    req = requests.get(url, proxies=proxies,
                       timeout=REQ_TIMEOUT, *args, **kwargs)
    req.raise_for_status()
    count = 0
    while self._meta_redirect(req.text):
      redirect_url = urljoin(url, self._meta_redirect(req.text))
      req = requests.get(redirect_url, proxies=proxies,
                         timeout=REQ_TIMEOUT, *args, **kwargs)
      req.raise_for_status()
      count += 1
      if count > 200:
        raise IOError("Too many redirects.")
    return req

  def _meta_redirect(self, content):
    '''Test whether a given soup object contains HTML redirect tag.'''
    soup   = BeautifulSoup(content, 'lxml')
    result = soup.find("meta", attrs={"http-equiv":re.compile("refresh", re.I)})
    if result:
      (wait, text) = result["content"].split(";")
      if text.lower().startswith("url="):
        url=text[4:]
        return url
    return None
  
  def _get_soup(self, params):

    '''Download run info according to the given params, parse and
    return a BeautifulSoup object.'''

    req = self._get_content(self.conf.sanger_lims_url,
                           params  = params)

    soup = BeautifulSoup(req.text, 'lxml')

    return(soup)

  def _consume_soup(self, soup):
    '''Parse the top-level NPG run listing and process each
    recently-changed row.'''
    tabs = soup.find_all('table')

    assert(len(tabs) == 1)

    tab = HtmlTable(tabs[0])

    # Find header.
    header = None
    for row in tab:
      elems = row.find_all('th')
      if len(elems) > 0:
        header = [ x.text for x in elems ]
        break

    assert(header is not None)
    try:
      datecol   = header.index('Status Date')
      statuscol = header.index('Status')
      instcol   = header.index('Instrument')
    except ValueError:
      raise ValueError('"Status", "Status Date" and/or "Instrument" columns not found in table header: '
                       % ", ".join(header))

    runnum_re = re.compile('_(\d+)$')

    # Parse rows.
    for row in tab:
      elems    = row.find_all('td')
      statdate = datetime.strptime(elems[datecol].text, DATEFORMAT)

      # Keep a record of the latest status date seen, and store that
      # in a cache file. This should help address time sync and timezone
      # issues.
      if statdate > self.running_status_date:
        self.running_status_date = statdate

      # If row status date is more recent than our last run, check to
      # see if it's of interest.
      if statdate > self.last_status_date:
        try:
          urlcol = header.index('Name')
        except ValueError:
          raise ValueError('"Name" column not found in header: ' % ", ".join(header))
        href    = elems[urlcol].find('a')
        url = urljoin(self.conf.sanger_lims_url, href['href'])
        try:
          (libdict, is_paired) = self._consume_recent_row(url)
        except requests.exceptions.HTTPError, err:
          LOGGER.warning("Unable to retrieve Run URL, skipping (%d error): %s",
                         int(err.response.status_code), url)
          continue

        # Pull out a few more details.
        statstr    = elems[statuscol].text
        instrument = elems[instcol].text
        runmatch = runnum_re.search(elems[urlcol].text)
        if runmatch:
          runnum = int(runmatch.group(1))
        else:
          raise ValueError("Unable to parse Run number from run name (%s)."
                           % elems[urlcol].text)
        facility = Facility.objects.get(code=self.conf.sanger_facility_code)
        flowcell = facility.code + 'run' + str(runnum)
        
        # Create or update libraries from libdict in the repository.
        for (flowlane, libcodes) in libdict.iteritems():
          for libcode in libcodes:

            # Libraries may not yet have been imported from our
            # inventory spreadsheet.
            try:
              library = Library.objects.get(code__iexact=libcode)
            except Library.DoesNotExist, err:
              LOGGER.warning("Library not available in repository: %s", libcode)
              self._missing_libcodes.add(libcode)
              continue

            # Update or create lane object.
            try:
              lane = Lane.objects.get(facility=facility,
                                      library=library,
                                      flowcell=flowcell,
                                      flowlane=flowlane)
              LOGGER.info("Updating lane for library %s, flowcell %s, lane %s",
                          library, flowcell, str(flowlane))
            except Lane.DoesNotExist, err:

              # Create new lane object.
              lane = Lane(facility=facility,
                          library=library,
                          flowcell=flowcell,
                          flowlane=flowlane,
                          lanenum=Lane.objects.next_lane_number(library))
              LOGGER.info("Creating new lane for library %s, flowcell %s, lane %s",
                          library, flowcell, str(flowlane))

            # Update lane object with other attrs.
            try: # Test that we have the authority to change repository lane status.
              auth = lane.status.authority == lane.facility
            except Status.DoesNotExist, err: # No status implies a new object.
              auth = True
            if auth:
              lane.status  = Status.objects.get(code=statstr,
                                                authority=lane.facility)
              lane.machine = Machine.objects.get(code__iexact=instrument)
              lane.rundate = statdate
              lane.paired  = is_paired

            # Write to database
            lane.save()

  def _study_is_ours(self, study_id, apiparts):
    study_url = "%s://%s/studies/%s.xml" % (apiparts.scheme, apiparts.netloc, str(study_id))
    req = self._get_content(study_url)
    root = ET.fromstring(req.text.encode('utf-8'))
    for owner in root.findall('./owners/owner'):
      login = owner.find('./login').text
      if login in self.conf.sanger_project_owners:
        return True
    return False

  def _consume_recent_row(self, url):
    '''Process a row from the main NPG run listing page.'''

    # Here we use a bit of a trick in which we access the AJAX URLs
    # addressed by the javascript embededd in the run webpage. The
    # alternative seemed to be to use e.g. PyQt4 and webkit to
    # actually run that javascript, but this seems much simpler and
    # probably not much more prone to breakage.
    req  = self._get_content(url)
    soup = BeautifulSoup(req.text, 'lxml')
    is_paired = self._is_paired(soup)
    libdict   = {}

    # Find the batch href and add a ".xml" extension. Without the
    # extension is password-protected; with the extension is freely
    # available. Go figure.
    #
    # We occasionally see failure here that appears to be due to
    # flakiness on the Sanger web server. We leave it for now and hope
    # that the next run is less flaky.
    batches = soup.find('th', text=re.compile('^Batch id', re.I))
    if batches is not None:
      batches = batches.next_sibling()
    if batches is None or len(batches) == 0:
      LOGGER.debug('Unable to parse batch URL from main NPG run listing. Skipping table row.')
      return (libdict, is_paired)
    batch = batches[0]
    batch_url = batch['href'] + '.xml'
    apiparts  = urlparse(batch_url)

    # Pull down the batch XML and extract library names, checking that
    # they belong to one of our studies.
    req = self._get_content(batch_url)
    root = ET.fromstring(req.text.encode('utf-8'))
    for lane in root.findall('./lanes/lane'):
      flowlane = int(lane.get('position'))
      libs = []
      for sample in lane.findall('./pool/sample'):
        study_id = sample.get('study_id')

        # We cache the study info to reduce the number of http queries.
        if study_id not in self._study_cache:
          self._study_cache[study_id] = self._study_is_ours(study_id, apiparts)

        if self._study_cache[study_id]:

          # Library is one of ours.
          lib_id  = sample.get('library_id')
          libnames = sample.get('library_name').split()
          assert(libnames[1] == lib_id)
          libs.append(libnames[0])
      libdict[flowlane] = libs

    return (libdict, is_paired)

  def _is_paired(self, soup):
    '''Pull out tags at the flowcell level.'''

    # PE/SE only at the moment.
    tagdiv = soup.find('div', attrs={ 'id':'current_tags' })
    if tagdiv is not None:
      return ('paired_read' in tagdiv.text.split())
    else:
      return False # A dubious default (FIXME review this).

  @transaction.atomic
  def _get_runs(self, querynum=100):

    # FIXME it might also be good to provide the option to pass in a
    # run ID to pull down specific information (e.g. for a MiSeq run
    # which would otherwise be ignored).
    params = { 'id_run_status_dict'   : 'all',
               'len'                  : querynum,
               'start'                : 0,
               'id_instrument_format' : 'all' }

    # Pull down an initial listing; we use no run data from this yet.
    soup = self._get_soup(params)

    # Figure out from the initial html what the appropriate
    # id_run_status_dict and id_instrument_format integers are.
    status_types = self._run_status_mapping(soup)
    instruments  = self._instrument_mapping(soup)

    facility   = Facility.objects.get(code=self.conf.sanger_facility_code)
    key_stages = Status.objects.filter(authority=facility).order_by('sortcode')

    # Cycle through queries for HiSeq, GAIIx (with an option on MiSeq)
    # for id_instrument_format.  These terms need to match the terms
    # as displayed on the main NPG runs listing webpage. If they don't
    # we will get a KeyError later on.
    for inst in self.conf.sanger_instruments:
      params['id_instrument_format'] = instruments[inst]
      LOGGER.info("Looking for %s runs...", inst)

      # As a sub-loop, cycle through id_run_status_dict for 'run
      # pending', 'run in progress', 'run complete', 'run archived',
      # 'qc complete'.
      for stage in key_stages:
        params['id_run_status_dict'] = status_types[stage.code]
        LOGGER.info("Querying Status: %s", stage.code)

        soup = self._get_soup(params)
        self._consume_soup(soup)

    if len(self._missing_libcodes) > 0:
      self._generate_missing_lib_email()

  def _generate_missing_lib_email(self):
    body  = "The following libraries being processed at Sanger"
    body += " have not yet been entered into our repository:\n\n"
    body += "\n".join(sorted(self._missing_libcodes))
    email_admins('[SangerPipe] Missing Libraries', body)

  def _generate_error_email(self, message):
    email_admins('[SangerPipe] Error Encountered', message)

  def get_runs(self, querynum=100):
    try:
      self._get_runs(querynum)
    except Exception, err:
      LOGGER.info("Rolled back database transaction.")
      self._generate_error_email(str(err))
      raise # must do this to signal failure to our context manager.
    
if __name__ == '__main__':

  with SangerLimsManager() as lims:
    lims.get_runs()
