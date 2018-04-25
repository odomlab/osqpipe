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

import os
import sys
from shutil import rmtree
from osqpipe.pipeline.laneqc import LaneFastQCReport
from osqutil.utilities import transfer_file, rezip_file, set_file_permissions
from osqutil.config import Config

from subprocess import Popen, PIPE

CONFIG = Config()

def run_qc(fnames, workdir, destination=None, cleanup=True, register=False):
    
    with LaneFastQCReport(fastqs=fnames, workdir=workdir, lane=0) as qc:
        # Generate qc reports
        qc.run_fastqc(qc.fastqs)
        qc.postprocess_results(qc.fastqs)

        # create list of disk files and if needed compress some of them before.
        dfiles = []
        # NB! This is not elegant, a better way of doing it would be if ftype.gzip and os.path.splitext(fname)[1] != CONFIG.gzsuffix:,
        #     However, this code is set up no to directly interact with database.
        for fn in qc.output_files:
            if fn.endswith('txt') or fn.endswith('tar'):
                dfn = rezip_file(fn)
                dfiles.append(dfn)
            else:
                dfiles.append(fn)

        if destination is not None:
            # transfer files to destination
            for dfn in dfiles:
                # set permissions
                set_file_permissions(CONFIG.group, dfn)
                # transfer file
                transfer_file(dfn, destination)
                
        if register:
            # register QC files in repository
            argslist = []
            for (fn,md5) in zip(qc.output_files, qc.output_md5s):
                argslist.append(os.path.basename(fn))
                argslist.append(md5) 
            # register files in repository
            cmd = "cs_addFile.py --qcfile -M --program_name %s " % qc.program_name
            cmd += " ".join(argslist)
            print "Executing \"%s\" ..." % cmd
            subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            (stdout, stderr) = subproc.communicate()
            retcode = subproc.wait()
            if stdout:
                sys.stdout.write(stdout)
            if stderr:
                sys.stderr.write(stderr)
            
        if cleanup:
            # remove local files
            # assuming fastqc report dir is still around, construct dirname.
            # NB! A cleaner way would be to save the dir name to self.bpath in postprocess_results in LaneQCReport class and use this value.
            #     Even better, perhaps LaneFastQCReport should be implemented to keep track of all temporary files it creates.
            for dfn in dfiles:
                os.remove(dfn)
                if dfn.endswith('pdf'):
                    fqc_dirname = os.path.splitext(dfn)[0]
                    rmtree(fqc_dirname)
                    zipfile = fqc_dirname + '.zip'
                    os.remove(zipfile)

if __name__ == '__main__':

    import argparse
      
    PARSER = argparse.ArgumentParser(
        description='Computes FastQC reports for a set of fastq files. Optionally, can move the files to destination and register in repository.')
    PARSER.add_argument('fnames', metavar='<files>', type=str, nargs='*',
                        help='List of files.')
    PARSER.add_argument('--workdir',dest='workdir', type=str,
                        help='Dir where the files are generated.')
    PARSER.add_argument('--cleanup',dest='cleanup', action='store_true', default=False,
                        help='Remove local files after being transferred to destination.')
    PARSER.add_argument('--destination',dest='destination', type=str, default=None,
                        help='Move files to destination (local folder or foreign destination user@host:/dir/)')
    PARSER.add_argument('--register',dest='register', action='store_true', default=False,
                        help='Register QC files in repository.')
    
    ARGS = PARSER.parse_args()

    if len(ARGS.fnames) < 1:
        PARSER.print_help()
        sys.exit(1)
        
    if ARGS.register and not ARGS.destination:
        sys.stderr.write("\n\nFiles destination not specified! Use --destination!\n\n")
        PARSER.print_help()
        sys.exit(1)

    if ARGS.cleanup and ARGS.destination is None:
        sys.stderr.write("Overriding --clean as --destination has not been specified\n")
        ARGS.cleanup = False
            
    run_qc(ARGS.fnames, workdir=ARGS.workdir, destination=ARGS.destination, cleanup=ARGS.cleanup, register=ARGS.register)
