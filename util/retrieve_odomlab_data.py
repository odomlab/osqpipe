#!/usr/bin/env python

'''
Simple client for Odom Lab repository REST API which can be used to
download data files to which the user has access.
'''

import requests
import json

def synchronise_datafiles(project=None, download=True, checksum=True,\
                            baseurl='http://localhost:8000/repository'):

  baseurl   = baseurl.rstrip('/')

  # First retrieve the user token, to be embedded in all future requests.
  tokenresp = requests.post("%s/api-token-auth/" % baseurl,
                            data={ 'username': user, 'password': password })
  assert(tokenresp.status_code == 200)
  token     = json.loads(tokenresp.content)['token']

  # Create a session object for future requests.
  sess = requests.Session()
  sess.headers.update({'Authorization': "Token %s" % token})

  # Now we retrieve some actual metadata.
  url    = "%s/api/projects" % baseurl
  resp   = sess.get(url)
  assert(resp.status_code == 200) ## FIXME detect token expiry here and refresh if necessary.
  output = json.loads(resp.content)

if __name__ == '__main__':

  from argparse import ArgumentParser

  parser = ArgumentParser(description=\
    'Download Odom lab data files with optional MD5 checksum confirmation.')

  parser.add_argument('--project', dest='project', type=str, required=False,
                      help='The optional name of the project to which downloads'
                      + ' will be restricted')

  parser.add_argument('--without-downlaod', dest='no_download', action='store_false',
                      help='Omit the file download step; this may be used to execute'
                      + ' checksums without repeating the file download.')

  parser.add_argument('--without-checksum', dest='no_checksum', action='store_false',
                      help='Omit the file MD5 checksum confirmation step. This is'
                      + ' provided as an option to skip what can be a disk I/O and'
                      + ' computationally intensive step.')

  parser.add_argument('--base-url', dest='baseurl', type=str,
                      default='https://dolab-srv003.cri.camres.org/django/repository',
                      help='The base repository URL to use for queries.')

  args = parser.parse_args()

  synchronise_datafiles(project  = args.project,
                        download = not args.no_download,
                        checksum = not args.no_checksum
                        baseurl  = args.baseurl)
