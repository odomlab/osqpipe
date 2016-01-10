#!/usr/bin/env python

'''
Simple client for Odom Lab repository REST API which can be used to
download data files to which the user has access.
'''

# NOTE that this script must contain *no dependencies* on
# osqpipe/osqutil modules.

import os
import requests
import json
import logging
import gzip
import hashlib
from contextlib import contextmanager

LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.handlers[0].setFormatter(\
  logging.Formatter("[%(asctime)s]dolab_%(levelname)s: %(message)s"))
LOGGER.setLevel(logging.DEBUG)

VERIFY_SSL_CERT=False

################################################################################
class ApiSession(requests.Session):
  '''
  Session subclass which will automatically renew expired
  authentication tokens.
  '''
  def __init__(self, username, password,
               base_url='http://localhost:8000/repository',
               *args, **kwargs):
    super(ApiSession, self).__init__(*args, **kwargs)
    self._rest_base_url      = base_url.rstrip('/')
    self._rest_username      = username
    self._rest_password      = password

  def _update_token(self):

    # Retrieve the user token, to be embedded in all future requests.
    LOGGER.debug("Refreshing authorization token.")
    api_token_url = "%s/api-token-auth/" % self._rest_base_url
    tokenresp = requests.post(api_token_url,
                              verify=VERIFY_SSL_CERT,
                              data={ 'username': self._rest_username,
                                     'password': self._rest_password, })
    assert tokenresp.status_code == 200
    api_token = json.loads(tokenresp.content)['token']
    self.headers.update({'Authorization': "Token %s" % api_token})

  def get(self, url, verify=VERIFY_SSL_CERT, *args, **kwargs):
    resp = super(ApiSession, self).get(url, verify=verify, *args, **kwargs)

    ## Detect token expiry here and refresh if necessary.
    if resp.status_code == 401:
      self._update_token()
      resp = super(ApiSession, self).get(url, *args, **kwargs)

    return resp

  def api_metadata(self, url):
    resp = self.get(url)
    assert resp.status_code == 200
    return json.loads(resp.content)

  def rest_download_file(self, url, local_filename):

    # The stream=True parameter keeps memory usage low.
    resp = self.get(url, stream=True)
    assert resp.status_code == 200
    with open(local_filename, 'wb') as outfh:
      for chunk in resp.iter_content(chunk_size=1024):
        if chunk: # filter out keep-alive new chunks
          outfh.write(chunk)
    return local_filename

################################################################################

@contextmanager
def flexi_open(filename, *args, **kwargs):
  '''
  Simple context manager function to seamlessly handle gzipped and
  uncompressed files.
  '''
  if os.path.splitext(filename)[1] == '.gz':
    handle = gzip.open(filename, *args, **kwargs)
  else:
    handle = open(filename, *args, **kwargs)

  yield handle

  handle.close()

def confirm_file_checksum(fname, checksum, blocksize=65536):

  with flexi_open(fname, 'rb') as fileobj:
    hasher = hashlib.md5()
    buf = fileobj.read(blocksize)
    while len(buf) > 0:
      hasher.update(buf)
      buf = fileobj.read(blocksize)

    if hasher.hexdigest() != checksum:
      LOGGER.warning("Local file checksum disagrees with repository: %s",
                     fname)

class OdomDataRetriever(object):

  __slots__ = ('session', 'with_download', 'with_checksum', '_base_url')

  def __init__(self, download=True, checksum=True,
               base_url='http://localhost:8000/repository' ):

    import getpass
    username = raw_input('Username: ')
    password = getpass.getpass('Password: ')

    # Create a session object for future requests.
    self.session = ApiSession(base_url = base_url,
                              username = username,
                              password = password )

    self.with_download = download
    self.with_checksum = checksum
    self._base_url     = base_url

  def process_liburl(self, url):

    libdict = self.session.api_metadata(url)
    LOGGER.info("Retrieved metadata for library %s", libdict['code'])
    for laneurl in libdict['lane_set']:
      self.process_laneurl(laneurl)

  def process_laneurl(self, url):

    lanedict = self.session.api_metadata(url)
    LOGGER.info("Retrieved metadata for lane %s (flowcell %s)",
                lanedict['flowlane'], lanedict['flowcell'])
    for filedict in lanedict['lanefile_set']:
      if filedict['filetype'] == 'fastq':
        self.process_lanefile(filedict)

  def process_lanefile(self, filedict):

    dl_fname = filedict['filename_on_disk']
    if os.path.exists(dl_fname):
      LOGGER.info("Skipping download of pre-existing file %s", dl_fname)
    else:
      if self.with_download:
        LOGGER.debug("Starting file download: %s", dl_fname)
        self.session.rest_download_file(filedict['download'], dl_fname)
        LOGGER.info("Downloaded file %s", dl_fname)
      else:
        LOGGER.debug("Skipping download as directed: %s", dl_fname)

    if not os.path.exists(dl_fname):
      LOGGER.warning("Downloaded file appears to be missing: %s", dl_fname)

    if self.with_checksum:
      confirm_file_checksum(dl_fname, filedict['checksum'])

  def synchronise_datafiles(self, project=None):

    # Now we retrieve some actual metadata.
    rootdict = self.session.api_metadata("%s/api/" % self._base_url)
    projlist = self.session.api_metadata(rootdict['projects'])

    if project is not None:
      projlist = [ proj for proj in projlist
                   if proj['code'].lower() == project.lower() ]

    liburls = sorted(list(set([ url for proj in projlist
                                for url in proj['libraries'] ])))

    for url in liburls:
      self.process_liburl(url)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description=\
    'Download Odom lab data files with optional MD5 checksum confirmation.')

  PARSER.add_argument('--project', dest='project', type=str, required=False,
                      help='The optional name of the project to which downloads'
                      + ' will be restricted')

  PARSER.add_argument('--without-download', dest='no_download', action='store_true',
                      help='Omit the file download step; this may be used to execute'
                      + ' checksums without repeating the file download.')

  PARSER.add_argument('--without-checksum', dest='no_checksum', action='store_true',
                      help='Omit the file MD5 checksum confirmation step. This is'
                      + ' provided as an option to skip what can be a disk I/O and'
                      + ' computationally intensive step.')

  PARSER.add_argument('--base-url', dest='baseurl', type=str,
                      default='https://dolab-srv003.cri.camres.org/django/repository',
                      help='The base repository URL to use for queries.')

  ARGS = PARSER.parse_args()

  RETRIEVER = OdomDataRetriever(download = not ARGS.no_download,
                                checksum = not ARGS.no_checksum,
                                base_url = ARGS.baseurl)
  RETRIEVER.synchronise_datafiles(project  = ARGS.project)
