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

import sys
import os

# set up logger
from osqutil.setup_logs import configure_logging
from logging import WARNING
LOGGER = configure_logging(level=WARNING)

# import config
from osqutil.config import Config

# For insertion of lane info:
import django
# from osqpipe.models import Lane, Status, Library, Facility, Machine, Adapter, ArchiveLocation
from osqpipe.models import Lane, Library, Lanefile, QCfile, LaneQC, Filetype
from osqutil.utilities import parse_repository_filename, sanitize_samplename

# set up config
DBCONF = Config()

django.setup()

class IntegrityCheck(object):
    
    def __init__(self, fn, merged_file=False, library=False, lane=False):
        # Find lane for the file
        self.fn = fn
        self.fn_base = os.path.basename(fn)
        self.lanes = None
        self.merged_file = merged_file
        self.library = library
        self.lane = lane

        self.find_lanes()

    def find_lanes(self):

        '''Find list of lane(s) associated with library, lane number and facility'''

        if not self.library:
            (code, facility, lanenum, pipeline) = parse_repository_filename(self.fn_base)
            # Merged files have only code in their prefix meaning failure by parse_repository_filename() above.
            if self.merged_file or self.lane:
                code = self.fn_base.split('_')[0]
        else:
            code = fn
        if code is None or code == '':
            LOGGER.error("Unable to extract code from filename %s." % self.fn_base)
            sys.exit(1)
        if self.merged_file or self.lane or self.library:
            self.lanes = Lane.objects.filter(library__code=code)
        else:
            self.lanes = Lane.objects.filter(library__code=code,
                                                lanenum=lanenum,
                                                facility__code=facility)
        if len(self.lanes) == 0:
            LOGGER.error("No lane associated with code \'%s\'." % code)
            sys.exit(1)

    def parse_bam_flagstat(self, fname):

        '''Parses file containing output from samtools flagstat.'''

        r1count = r2count = rcount = 0
        properly_paired_percent = mapped_percent = ""

        if not os.path.exists(fname):
            LOGGER.error("%s missing!\n" % fname)
            return (r1count, r2count, rcount, mapped_percent, properly_paired_percent)

        fh = open(fname, 'rb')
        for line in fh:
            line = line.rstrip('\n')
            if line.endswith('read1'):
                r1count = int(line.split(' + ')[0])
                # print "Read1='%d'" % r1count
            if line.endswith('read2'):
                r2count = int(line.split(' + ')[0])
                # print "Read2='%d'" % r2count
            if line.endswith('paired in sequencing'):
                rcount = int(line.split(' + ')[0])
            if 'mapped (' in line:
                mapped_percent = line.split('(')[1].split('%')[0]
            if 'properly paired' in line:
                properly_paired_percent = line.split('(')[1].split('%')[0]
                break

        if (r1count + r2count) != rcount:
            SystemExit("Number of left (%d) and right (%d) reads not equal!\n" % (r1count, r2count))

        return (r1count, r2count, rcount, mapped_percent, properly_paired_percent)

    def check_bam_integrity(self, merged_bam=False):

        '''Compares number of reads in bam with number of reads in associated lane.
        Expects precomputed fn.bam.flagstat in path of fn.'''

        nr_of_reads = 0
        for lane in self.lanes:
            # print "Bam=%s\tLane:" % self.fn
            # print lane
            nr_of_reads += lane.passedpf
            if not merged_bam:
                if len(self.lanes) > 1:
                    LOGGER.error("Found multiple lanes for '%s': %s",
                                 self.fn, ", ".join([x.id for x in self.lanes]))
                break

        # Parse flagstat
        fn_flagstat = self.fn + '.flagstat'
        (r1count, r2count, rcount, mapped_percent, properly_paired_percent) = self.parse_bam_flagstat(fn_flagstat)

        # Compare number of reads with content of flagstat info for the bam.
        if r1count != nr_of_reads:
            LOGGER.error("Reads in bam=%d and database=%d not same!\t%s" % (r1count, nr_of_reads, self.fn) )
            print "File DIFFERENT!\t%s" % (fn)
            return "DIFFERENT"
        else:
            print "File OK!\t%s" % (fn)
            return "OK"

    def check_file_on_disk(self, fobj):
        '''Checks if file associated with file object is accessible on disk.'''

        root_path = DBCONF.repositorydir
        if fobj.archive:
            root_path = fobj.archive.root_path        
        fsuffix = ''
        if fobj.filetype.gzip:
            fsuffix = '.gz'
        fpath = os.path.join(root_path, fobj.libcode, fobj.filename + fsuffix)
        return os.path.isfile(fpath)

    def check_lane_fq_integrity(self, fobj):
        if not self.check_file_on_disk(fobj.filename):
            print "ERROR: %s FASTQ FILE MISSING on disk!" % fobj.filename
            if lanefileobj.filetype.code=='fq':
                expected_qcfiles[lanefile.filename + "_fastqc.pdf"] = 1
                expected_qcfiles[lanefile.filename + "_fastqc.tar"] = 1
                expected_qcfiles[lanefile.filename + "_fastqc.txt"] = 1
                
    def check_lane_integrity(self):
        
        '''Check integrity of lane associated records.'''
        for lane in self.lanes:
            nr_of_fqfiles = 0
            expected_qcfiles = dict()
            # check if recorded lanefiles exist on disk
            for lanefile in Lanefile.objects.filter(lane__id=lane.id):
                if not self.check_file_on_disk(lanefile):
                    print "ERROR: %s FASTQ FILE MISSING on disk!" % lanefile.filename
                if lanefile.filetype.code=='fq':
                    nr_of_fqfiles += 1
                    expected_qcfiles[lanefile.filename + "_fastqc.pdf"] = 1
                    expected_qcfiles[lanefile.filename + "_fastqc.tar"] = 1
                    expected_qcfiles[lanefile.filename + "_fastqc.txt"] = 1
            if lane.paired and nr_of_fqfiles!= 2:
                print "ERROR: Paired-end lane. Two files expected, %d found!" % nr_of_fqfiles
            if not lane.paired and nr_of_fqfiles!=1:
                print "ERROR: Single-end lane. One file expected, %d found!" % nr_of_fqfiles

            # check if recorded QC files exist on disk
            laneqcs = LaneQC.objects.filter(lane = lane)
            for lqc in laneqcs:
                qcfs = QCfile.objects.filter(laneqc=lqc)
                for qcf in qcfs:
                    if qcf.filename in expected_qcfiles:
                        del expected_qcfiles[qcf.filename]
                        if not self.check_file_on_disk(qcf):
                            print "ERROR: %s QCFILE MISSING on disk!" % qcf.filename
            # check if database records for any of the expected QC files were missing.
            for fname in expected_qcfiles:
                print "ERROR: %s MISSING in repository!" % fname

if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Check file integrity')

  PARSER.add_argument('fn', metavar='<filename>', type=str,
                      help='Filename.', nargs='+')

  PARSER.add_argument('-B', '--merged_bam', dest='merged_bam', action='store_true',
                      help='Merged bam file(s) across all lanes for the library. Expects <filename>.flagstat in the same path. Checks if number of reads in file(s) equals to reads for the library.')

  PARSER.add_argument('-b', '--bam', dest='bam', action='store_true',
                      help='Bam file(s). Expects <filename>.flagstat in the same path. Checks if number of reads in file(s) is the same as for lane in repository.')

#  PARSER.add_argument('-c', '--code', dest='code', action='store_true',
#                      help='Library code.')

  PARSER.add_argument('-l', '--lane', dest='lane', action='store_true',
                      help='A filename following Odom file name convention. Checks presence of lane associated files and records.')

  PARSER.add_argument('-L', '--library', dest='library', action='store_true',
                      help='A filename following Odom file name convention. Checks presence of library associated files and records.')  

  PARSER.add_argument('-a', '--alignment', dest='alignment', action='store_true',
                      help='A filename following Odom file name convention. Checks presence of library associated files and records.')  

  ARGS = PARSER.parse_args()

  for fn in ARGS.fn:
    icheck = IntegrityCheck(fn, merged_file=ARGS.merged_bam, library=ARGS.library, lane=ARGS.lane)
    
    if ARGS.merged_bam:
        icheck.check_bam_integrity(merged_bam=True)
    if ARGS.bam:
        icheck.check_bam_integrity()
    if ARGS.library:            
        icheck.check_lane_integrity()
        # icheck.check_library_integrity()
    if ARGS.lane:
        icheck.check_lane_integrity()
    if ARGS.alignment:
        sys.stderr.write("Integrity check for alignment not implemented.")
        # icheck.check_alignment_integrity()
