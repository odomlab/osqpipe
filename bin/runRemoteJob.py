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

'''Simple script to take a list of files and a command to run on them,
copy the files over to the currently configured remote job host
(either the cluster or a suitable alternative alignment desktop host),
and run the command. No attempt is made to retrieve the output; if
desired, such behaviour needs to be embedded in the remote command
itself.'''

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

from osqutil.config import Config
from osqpipe.pipeline.bwa_runner import DesktopJobSubmitter, ClusterJobSubmitter

def run_job(cmd, files, append=False, mem=2000, testmode=False):

  if files is None:
    files = []

  config = Config()

  try:
    host = config.althost
    assert(host != '')
    runner = DesktopJobSubmitter(test_mode=testmode)
  except Exception:
    runner = ClusterJobSubmitter(test_mode=testmode)

  if append:
    cmd = " ".join([cmd] + files)

  LOGGER.info("Transferring data files...")
  runner.transfer_data(files)

  LOGGER.info("Running command...")
  runner.submit_command(cmd, mem=mem)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description='Run a job on a remote host.')

  PARSER.add_argument('-c', '--command', dest='command', type=str, required=True,
                      help='The command to run on the remote host; if the \'-a\' option is'
                      + ' used then the filenames given here will be appended to this command'
                      + ' prior to execution.')

  PARSER.add_argument('-a', '--append', dest='append', action='store_true',
                      help='Append the filenames to the command. If this option is omitted '
                      + 'then it is assumed that the given command is sufficent.')

  PARSER.add_argument('files', metavar='<files>', type=str, nargs='+',
                      help='The name of the file or files to transfer for processing.')

  PARSER.add_argument('-m', '--memory', dest='memory', type=int, required=True, default=2000,
                      help='The amount of memory resource to request (if supported) on the'
                      + ' remote host. Used for example in submitting commands to LSF queues.')

  PARSER.add_argument('-t', '--test', dest='testmode', action='store_true',
                      help='Turn on test mode.')

  ARGS = PARSER.parse_args()
  run_job(ARGS.command, ARGS.files, ARGS.append, ARGS.memory, ARGS.testmode)

