#!/usr/bin/env python
#
# $Id$
#
# Script to:
# (a) parse the inventory spreadsheet
# (b) check for libraries not represented in the repository
# (c) insert new libraries into the repository
# (d) optionally update all metadata for libraries already found in the repository.

import os
import re
from xlrd import open_workbook

from osqutil.setup_logs import configure_logging
from logging import INFO, DEBUG
LOGGER = configure_logging(level=INFO)

# New in Django 1.7 and above.
import django
django.setup()

from osqutil.config import Config
from osqpipe.models import Library

from osqpipe.pipeline.library import LibraryHandler

class InventoryImporter(object):

  __slots__ = ('verbose', 'conf', 'libhandler', 'test_mode')

  def __init__(self, verbose=False, test_mode=False):
    self.verbose = verbose
    self.test_mode = test_mode
    self.conf    = Config()
    self.libhandler = LibraryHandler(interactive=False, fuzzy=True,
                                     test_mode=test_mode)
    if verbose:
      LOGGER.setLevel(DEBUG)
    else:
      LOGGER.setLevel(INFO)

  def throw_exception(self, message, exceptionType=Exception):
    LOGGER.error(message)
    raise(exceptionType(message))

  def get_latest_spreadsheet(self, path):

    '''Given a path, identify the newest .xls or .xlsx file in that
    directory (without recursion) and return the full path.'''

    LOGGER.info("Identifying latest ChIP_Inventory file ...")
    xlspat = re.compile(r'\.xlsx?$')
    latestFile = None
    latestTime = 0
    for item in os.listdir(path):
      if item[0] == '.': # skip dotfiles
        continue
      mo = xlspat.search(item)
      item = os.path.join(path, item)
      if os.path.isdir(item) or mo is None:
        continue
      t = os.path.getmtime(item)
      if t > latestTime:
        latestFile = item
        latestTime = t
    if latestFile is not None:
      latestFile = os.path.realpath(latestFile)
      LOGGER.info("File found: %s" % (latestFile,))
    else:
      self.throw_exception("No suitable file found.")
    return latestFile

  def get_latest_work_sheet(self, path, sheet):

    '''Given a path and a sheet name find the newest Excel file and
    look for the named sheet therein, returning the parsed object.'''
      
    xls = self.get_latest_spreadsheet(path)
    wb  = open_workbook(xls)
    for s in wb.sheets():
      if s.name == sheet:
        return s
    self.throw_exception("Unable to identify target sheet (%s) in workbook." % (sheet,))
    return None

  def get_work_sheet_iterator(self, path, sheet):

    '''Create an iterator object which can be used to walk through the
    list of libraries.'''

    shobj = self.get_latest_work_sheet(path, sheet)
    for rownum in range(shobj.nrows):
      yield shobj.row(rownum)

  def import_work_sheet(self, path, sheet='DO_list', minlib=None):

    row_iter = self.get_work_sheet_iterator(path, sheet)
    header = []

    # Detect the header row.
    space_re = re.compile(' ')
    for row in row_iter:
      rowvals = [ space_re.sub('', unicode(x.value).lower().strip()) for x in row ]
      if 'libraryid' in rowvals:
        header = rowvals
        break
    if len(header) == 0:
      raise ValueError("Header row not found in spreadsheet"
                       + " (check for LibraryID column).")

    # Very quick sanity check against our current spreadsheet file format.
    known = ('sort', 'libraryid', 'assaytype', 'experiment', 'genome',
             'strain', 'tissue', 'cellline', 'tissueprocessing', 'projects',
             'individual', 'sex', 'factor', 'antibody', 'lotnumber', 'condition',
             'barcode', 'barcode2', 'barcode(column)', 'barcode2(row)',
             'barcodetype', 'linkerset', 'pairedend', 'protocol',
             'notes', 'status', 'replicate', 'biopsyid', 'solexaid',
             'sequencingrunsatcri:', 'sequencingrunsatsanger:',
             'total#ofmappedreadsobtained:','comment',
             'mappedtogenomeversion:(e.g.hg18mm9usw)')
    for colname in header:
      if colname not in known:
        LOGGER.warning("Unrecognised column in input: %s", colname)

    # If desired, skip over most of the file without running expensive database queries.
    if minlib is not None:
      codere = re.compile('(\w+)(\d+)')
      minmo = codere.match(minlib)
      if minmo is None:
        raise ValueError("Starting library code does not fit our standard format (\w+\d+).")
      for row in row_iter:
        code = row[header.index('libraryid')]
        libmo = codere.match(code)
        if libmo is None:
          raise ValueError("Library code in sheet does not fit our standard format (\w+\d+).")
        if minmo.group(1) == libmo.group(1) and minmo.group(2) >= libmo.group(2):
          self.import_data_row([x.value for x in row], header)
          break

    # Read in the rest of the file.
    for row in row_iter:
      self.import_data_row([x.value for x in row], header)

  @staticmethod
  def extract_barcode_int(key, rowdict):
    int_re  = re.compile('(\d+)')
    barcode = None
    if key in rowdict:
      match = int_re.search(unicode(rowdict[key]))
      if match:
        barcode = match.group(1)
    return barcode

  @staticmethod
  def process_optional_values(rowdict):

    # Keys are normalised spreadsheet column names; values are
    # database column names (or at least, codes which can be
    # interpreted by the LibraryHandler class)..
    mapping = {'linkerset'  : 'linkerset',
               'experiment' : 'chipsample',
               'strain'     : 'strain',
               'individual' : 'individual',
               'factor'     : 'factor',
               'antibody'   : 'antibody',
               'pairedend'  : 'paired',
               'lotnumber'  : 'lot_number',
               'sex'        : 'sex',
               'comment'    : 'comment',
               }

    # Read in any optional data into optvals.
    optvals = {}
    for sscol in mapping:
      if sscol in rowdict and rowdict[sscol] != '':
        optvals[ mapping[sscol] ] = rowdict[sscol]

    return optvals

  def _munge_truseq_neb_barcode_info(self, rowdict, barcode):
    '''
    Specialised method to handle some of the niceties surrounding
    TruSeq adapter parsing. This method is designed to handle rows for
    which the protocol contains \b(truseq|neb)\b.
    '''
    adapter = linkerset = None

    prottag = rowdict['protocol'].lower()

    # TruSeq adapters have standard names in the repository.
    if rowdict['assaytype'].lower() in ('smrnaseq', 'rip-smrnaseq'):
      adapter = 'smRNA_TS_' + barcode   # smRNAseq

      # Often, linkerset is not given but it's easily guessed.
      if 'linkerset' not in rowdict or rowdict['linkerset'] == '':
        linkerset = 'TruSeqSmRNAIndex' + barcode

    elif ( int(barcode) > 500 and int(barcode) < 509 )\
         or ( int(barcode) > 700 and int(barcode) < 713 ):
      adapter = 'TS_D' + barcode   # dual indexing adapter set.

    else:
      adapter = 'TS_' + barcode    # not smRNAseq

      if int(barcode) > 12 and rowdict['protocol'].lower() == 'neb':
        LOGGER.error('NEB barcode index greater than 12 used; confirm sequences match TruSeq!')
        raise ValueError()

    return (adapter, linkerset)

  def munge_barcode_info(self, rowdict, libcode, code_column='barcode', optional=False):

    adapter = linkerset = None

    needs_adapter = False
    if 'barcodetype' in rowdict:
      needs_adapter = str(rowdict['barcodetype']).strip() != ''

    barcode = self.extract_barcode_int(code_column, rowdict)
    if code_column in rowdict:
      needs_adapter |= str(rowdict[code_column]).strip() != ''

    # Handle failures, inconsistencies and other corner cases.
    if barcode:  # If we have a barcode number we need to know the scheme

      # NEB and TruSeq barcodes are identical, at least up to TS_12.
      if 'protocol' in rowdict:
        prottag = rowdict['protocol'].lower()

        if re.search(r'\b(?:truseq ?lt)\b', prottag):
          adapter = 'TSLT_%d' % int(barcode)

        elif re.search(r'\b(?:truseq|neb)\b', prottag):
          (adapter, linkerset) = self._munge_truseq_neb_barcode_info(rowdict, barcode)

        elif re.search(r'\b(?:nextera|nextera ?xt)\b', prottag):
          adapter = 'NXT_N' + barcode

        elif re.search(r'\b(?:thruplex)\b', prottag):
          if int(barcode) < 150:
            adapter = 'iPCRtagT' + barcode
          elif ( int(barcode) > 500 and int(barcode) < 509 )\
               or ( int(barcode) > 700 and int(barcode) < 713 ):
            adapter = 'TP_D' + barcode   # dual indexing adapter set, as for TruSeq.
          else:
            LOGGER.error('Uncertain which ThruPlex adapters have been used: %s',
                         libcode)
            raise ValueError()

        elif re.search(r'\b(?:haloplex|agilent ?haloplex)\b', prottag):
          adapter = 'HAL' + barcode

        # FIXME the following protocol list should be simplified to
        # reflect whatever we end up actually using.
        elif re.search(r'\b(?:agilent ?sureselect ?xt|sureselect ?xt|sureselect)\b', prottag):

          if int(barcode) > 16:
            LOGGER.error('SureSelectXT barcode index greater than 16 used.')
            raise ValueError("Unlikely SureSelectXT barcode used: %s" % barcode)

          # This protocol also supports the use of 8bp adapters named
          # A01 through H12 (for 96-well plates). Integers are not
          # going to be sufficient to disambiguate.
          if re.match('[A-H]', rowdict[code_column]):
            barcode = rowdict[code_column]

          adapter = 'XT_' + barcode

        elif 'barcodetype' in rowdict:
          # In some cases we can fall back to barcodetype
          bctag = rowdict['barcodetype'].lower()
          if re.search(r'\b(?:truseq|neb)\b', bctag):
            (adapter, linkerset) = self._munge_truseq_neb_barcode_info(rowdict, barcode)
          else:
            LOGGER.error('Uncertain which adapter scheme (e.g. TruSeq) has been used: %s',
                         libcode)
            raise ValueError()
          
      else:
        LOGGER.error('Protocol not specified in spreadsheet for library %s', libcode)
        raise ValueError()
      
    elif needs_adapter and not optional:

      # If there's a non-empty barcode we've failed to parse, that's a problem.
      LOGGER.error('Unable to guess barcode from %s/barcodetype columns.', code_column)
      raise ValueError()

    return (adapter, linkerset)

  def process_project_codes(self, rowdict):

    projcodes = []
    projcol = 'projects'
    if projcol in rowdict and rowdict[projcol] != '':
      # The regex in the following line defines the list of valid
      # characters in a project code. Note that whitespace and
      # forward-slash will break the website url construction, so
      # avoid them.
      projcodes = [ x.strip() for x in re.split('\W+', rowdict[projcol]) ]
    return projcodes
  
  def import_data_row(self, row, header):
    '''Given a list of row data and a header list, attempt to load it
    into the database.'''

    rowdict = dict(zip(header, row))
    if 'libraryid' not in rowdict or rowdict['libraryid'] == '':
      LOGGER.warning("Skipping row with no Library ID code.")
      return
    # We've standardized on lowercase library codes. Note that for
    # this to work with formula-generated libcodes, the libraryId
    # column must be specified as text format.
    libcode = rowdict['libraryid'].lower()
    LOGGER.debug("Checking for library ID %s in database...", libcode)

    # Actually handle the data here.
    if Library.objects.filter(code__iexact=libcode).count() > 0:
      LOGGER.debug("Found library ID in database already; skipping: %s", libcode)
      return  # FIXME we may wish to update as well.

    # Check our required arguments
    for key in ('assaytype', 'genome'):
      if key not in rowdict or rowdict[key] == '':
        LOGGER.error("Library does not have required %s annotation: %s",
                     key, libcode)
        return

    if len([ key for key in ('tissue','cellline')
             if key not in rowdict or rowdict[key] == '' ]) < 1:
      LOGGER.error("Either tissue or cell line must be specified: %s", libcode)
      return

    optvals = self.process_optional_values(rowdict)

    projcodes = self.process_project_codes(rowdict)

    # Overwrite tissue with cellline if the latter is given.
    if 'cellline' in rowdict and rowdict['cellline'] != '':
      tissue = rowdict['cellline']
    else:
      tissue = rowdict['tissue']

    # Condition is occcasionally supplied; if so, we record it.
    if 'condition' in rowdict and rowdict['condition'] != '':
      optvals['condition'] = rowdict['condition']

    # Munge the barcode/barcodetype/protocol info into an adapter string.
    try: # Test New vs. Old column naming, in case we ever switch back.

      (optvals['adapter'],  optvals['linkerset']) = \
          self.munge_barcode_info(rowdict, libcode, code_column='barcode(column)')
      if optvals['adapter'] is None:
        (optvals['adapter'],  optvals['linkerset']) = \
            self.munge_barcode_info(rowdict, libcode, code_column='barcode')

      (optvals['adapter2'], discarded) = \
          self.munge_barcode_info(rowdict, libcode, code_column='barcode2(row)', optional=True)
      if optvals['adapter2'] is None:
        (optvals['adapter2'], discarded) = \
            self.munge_barcode_info(rowdict, libcode, code_column='barcode2', optional=True)
      
      if discarded is not None:
        raise ValueError("Unexpectedly identified a linkerset (%s) while parsing adapter2." % discarded)

    # Here we finally admit defeat.
    except ValueError, err:
      LOGGER.error("Barcode parsing failed for library %s. Skipping.", libcode)
      return

    # We handle fuzzy matching in the LibraryHandler class.
    self.libhandler.add(libtype = rowdict['assaytype'],
                        code    = libcode,
                        genome  = rowdict['genome'],
                        tissue  = tissue,
                        projcodes = projcodes,
                        opts    = optvals)


if __name__ == '__main__':

  import argparse
  from logging import FileHandler, Formatter

  PARSER = argparse.ArgumentParser(
    description='Update the repository with information'
    + ' from the inventory spreadsheet.')

  PARSER.add_argument('-d', '--dir', dest='dir', type=str, required=True,
                      help='Directory containing the inventory spreadsheets.')

  PARSER.add_argument('-s', '--sheet', dest='sheet', type=str, default='DO_list',
                      help='The worksheet name (default: DO_list).')

  PARSER.add_argument('-t', '--test', dest='test_mode', action='store_true',
                      help='Test mode (no changes to the database).')

  PARSER.add_argument('-l', '--logfile', dest='logfile', type=str, required=False,
                      help='(Optional) Log file name to store output'
                      + ' (appends to an existing file).')

  PARSER.add_argument('--start-from', dest='minlib', type=str, required=False,
                      help='(Optional) Library code from which to commence import (to save time).')

  ARGS = PARSER.parse_args()

  if ARGS.logfile:
    HND  = FileHandler(ARGS.logfile)
    FRMT = Formatter("[%(asctime)s]%(levelname)s: %(message)s")
    HND.setFormatter(FRMT)
    LOGGER.addHandler(HND)

  IMP = InventoryImporter(test_mode=ARGS.test_mode)

  IMP.import_work_sheet(path=os.path.realpath(ARGS.dir), sheet=ARGS.sheet, minlib=ARGS.minlib)
