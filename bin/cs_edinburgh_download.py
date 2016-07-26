#!/usr/local/bin/python
#
# $Id$

""" Run Data Download from Edinburgh """
__author__ = "Margus Lukk"
__date__ = "13 July 2016"
__version__ = "0.2"

import sys
import os
import re
from subprocess import Popen, PIPE

# set up logger
from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

# import config
from osqutil.config import Config

# set up cluster job submitter
from osqutil.cluster import ClusterJobSubmitter

# Import some functions from utilities package
from osqutil.utilities import call_subprocess, _checksum_fileobj # Note that for _checksum_fileobj this may not work as _ indicates that the function is not part of API!

# For insertion of lane info:
import django
from osqpipe.models import Lane, Status, Library, Facility, Machine, Adapter, ArchiveLocation

# set up config
DBCONF = Config()

## In the heart of the download code lies following:
# 1. Given [USER ID], construct links to all files for this [userid]. NB! What edinburgh calls userid will be in our case donumber.
# 2. Download file and download md5sum
# 3. Compare file and md5sum, if not same, go back to 2. Keep track how many attempts file was downloaded (limit to 10.)
# 4. Set file on diks .done
# 5. Get lane info for lane and flowcell. If none, register flowcell.
# 6. Set flowcell status 'downloaded'
## There are multiple entry points to this code

## Note that the code depends on following variables defined in config:
# acredfile - containing credentials for connecting to Edinburgh
# ahost - host of the aspera server in Ediburng (edgen-dt.rdf.ac.uk)
# apport - aspera port for TCP communication
# aoport - aspera port for UDP communication
# arate - aspera download rate. E.g. 500M = 500Mbit/s
# athreads - number of parallel downloads to execute

## TODO:
# 1. Test the script with -l and -i options. The script has been tested with -f and -F options.
# 2. Test the script with -i option with different number of threads. See the athread option in the config file.
# 3. Check that the logging of failed download commands to a file specified in config as 'faileddownloads' is working fine.

def read_credentials(credentials_file):
    '''Reads credentials from a file to memory.'''
    # Assumes following file structure
    # username\tuser str
    # password\tpassword str
    # 
    credentials = {'username':None, 'password':None}
    fh = open(os.path.join(os.environ['HOME'], credentials_file),'rb')
    for line in fh:
        line = line.rstrip('\n')
        cols = line.split('\t')
        if cols[0] in credentials:
            credentials[cols[0]] = cols[1]
    for key in credentials:
        if credentials is None:
            sys.exit('Ill formated credentials file. No value for \'%s\'!\n\n' % key)
    return credentials

def compute_md5(fn):

    md5 = None

    LOGGER.info("Computing md5 for %s ...", fn)
    with open(fn, 'rb') as fileobj:
        md5 = _checksum_fileobj(fileobj)
    return md5

def run_command(cmd, shell=False, path=None):
    '''A simple wrapper for running a command cmd for cases where using utilities.call_subprocess is not appropriate.
    Agnostic whether running command was successful. Returns list of length 3 containing stdout, stderr and retcode'''
    
    # Set our PATH environmental var to point to the desired location.                                                                                          
    oldpath = os.environ['PATH']
    if path is not None:
        if type(path) is list:
            path = ":".join(path)
        os.environ['PATH'] = path
    else:
        LOGGER.warn("Subprocess calling external executable using undefined $PATH.")

    subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=shell)
    (stdout, stderr) = subproc.communicate()
    retcode = subproc.wait()

    # switch environment back to what it was before.
    os.environ['PATH'] = oldpath

    return (stdout, stderr, retcode)

def parse_read_header(fn):
    '''Extracts flowcell from the name of the first read in the file. Assumes either uncompressed, gzipped or bzipped fastq input file'''
    
    # Construct command for extraction of flowcell from file
    cmd = 'cat %s' % fn
    if fn.endswith('.gz'):
        cmd = 'zcat %s' % fn
    if fn.endswith('.bz2'):
        cmd = 'bzcat %s' % fn
    cmd = cmd + ' | head -n 1'
    
    # Run command
    (stdout, stderr, retcode) = run_command(cmd, shell=True)
    
    # stderr is likely to contain error trown by zcat or bzcat due to head existing after reading only one line.
    # e.g. 'gzip: stdout: Broken pipe' or 'bzcat: I/O or other error, bailing out. Possible reason follows.'
    # Assume Illumina default read header:
    m = re.search('^@([0-9A-Za-z]+):\d+:([0-9A-Za-z]+):(\d):.*$', stdout)
    if m:
        machine = m.group(1)
        flowcell = m.group(2)
        flowlane = m.group(3)
        if len(flowcell) >= 7:
            return (machine,flowcell,flowlane)
    sys.exit("Failed extracting flowcell information from file '%s'\n" % fn)

class edFile(object):
    def __init__(self, edstem, fname, md5, pair):

        # Stem of the edinburgh filename
        # E.g for '2016-07-12/X15060P001B01/raw_data/160704_E00328_0103_AHNCMKCCXX_4_ATTACTCG_R2.fastq.gz'
        # The stem would be '2016-07-12/X15060P001B01/raw_data/160704_E00328_0103_AHNCMKCCXX_4_ATTACTCG'
        self.edstem = edstem
        
        # File names in edinburgh aspera server with path.
        # E.g. '2016-07-12/X15060P001B01/raw_data/160704_E00328_0103_AHNCMKCCXX_4_ATTACTCG_R2.fastq.gz'
        self.file1 = None
        self.file2 = None
        self.file1_md5 = None 
        self.file2_md5 = None

        # id and status of the associated lane
        self.laneid = None
        self.status = 'unknown'
        
        # File name stem for this lane in our repository
        self.repstem = None
        # File names for pair1 and pair2
        self.rep_file1 = None
        self.rep_file2 = None
            
        self.addFile(fname, md5, pair)

    def addFile(self, fname, md5, pair):

        if pair == 1:
            if self.file1 is not None:
                LOGGER.error("%s already recorded! Duplicate entry? Exiting!")
                sys.exit(1)
            else:
                self.file1 = fname
                self.file1_md5 = md5
        elif pair == 2:
            if self.file2 is not None:
                LOGGER.error("%s already recorded! Duplicate entry? Exiting!")
                sys.exit(1)
            else:
                self.file2 = fname
                self.file2_md5 = md5
        else:
            LOGGER.error("File pair can not be %s. Either 1 or 2 expected!", pair)
            sys.exit(1)


class ed_data_handler(object):
    def __init__(self, project, maxattempts=5, threads=None, aname='ebiark'):

        # Set the Edinburgh aspera project directory.
        self.project = project

        # Set maximum number of attempts to download any single file
        self.maxattempts = maxattempts
    
        # Set the name of the archive where the fastq files should be saved
        self.aname = aname
    
        # Get following values from conf file
        self.ahost = DBCONF.ahost
        self.aPport = DBCONF.apport
        self.aOport = DBCONF.aoport
        self.arate = DBCONF.arate    
        if threads is not None:
            self.athreads = threads
        else:
            self.athreads = int(DBCONF.athreads)
        self.path = DBCONF.clusterpath

        # Create project directory in default incoming directory
        self.destination = os.path.join(DBCONF.incoming, self.project)
        if not os.path.isdir(self.destination):
            try:
                os.makedirs(self.destination)
            except OSError:
                LOGGER.error("Failed to create project directory %s", self.destination)
                sys.exit(1)
        LOGGER.info("Setting project directory to %s", self.destination)
        self.failedcommands = os.path.join(self.destination, DBCONF.faileddownloads)

        # Load aspera credentials
        self.credentials = read_credentials(DBCONF.acredfile)
        os.environ['ASPERA_SCP_PASS'] = self.credentials['password']        

        # Define project files and containers holding project related info
        self.sample2library = dict()        
        self.md5sums_fn = 'all_md5sums.txt'
        self.summary_fn = 'summary_metrics.csv'
        self.report_fn = 'project_%s_report.pdf' % self.project
        # Define lists to hold tri information about each edinburgh file: md5sum, filename and process status
        self.edfiles = dict()

    def _ed_parse_project_files(self):
        # If not yet downloaded, download files in project root
        pfiles = [self.md5sums_fn, self.report_fn, self.summary_fn]
        for pfile in pfiles:
            res = edd.ed_get_file_with_md5(pfile, None, download_md5=False)
            if res != 0:
                LOGGER.error("Unable to download %s. Exiting!", pfile)
                sys.exit(1)

        # Parse info about libraries and samples
        with open(os.path.join(self.destination,self.summary_fn), 'rb') as fh:
            firstline = True
            for line in fh:
                if firstline:
                    firstline = False
                    continue
                line = line.rstrip('\n')
                cols = line.split()
                self.sample2library[cols[1]] = cols[2].lower()
        fh.close()

        # Parse info from all_md5sums file
        with open(os.path.join(self.destination,self.md5sums_fn), 'rb') as fh:
            for line in fh:
                line = line.rstrip('\n')
                (md5, fname) = line.split()
                if fname.endswith('_R1.fastq.gz'):
                    edstem = fname.rstrip('_R1.fastq.gz')
                    pair = 1
                elif fname.endswith('_R2.fastq.gz'):
                    edstem = fname.rstrip('_R2.fastq.gz')
                    pair = 2
                else:
                    LOGGER.error("%s suffix not recognised! '_R?.fastq.gz' expected!")
                    sys.exit(1)
                if edstem not in self.edfiles:
                    self.edfiles[edstem] = edFile(edstem, fname, md5, pair)
                else:
                    self.edfiles[edstem].addFile(fname, md5, pair)
        fh.close()
        
        # Sanity check that all samples/libraries have files for both pairs.
        for edstem in self.edfiles:
            if self.edfiles[edstem].file1 is not None and self.edfiles[edstem].file2 is not None:
                continue
            else:
                LOGGER.error("File missing in pair \'%s\' and \'%s\'", self.edfiles[edstem].file1, self.edfiles[edstem].file2)
                sys.exit(1)
                        
    def create_lane(self, lib, flowcell, flowlane, facility, status, rundate, sampleid, runnumber, instrument, notes):

        try:
            lane = Lane.objects.get(flowcell=flowcell, flowlane=flowlane, library=lib)    
        except Lane.DoesNotExist, _err:            
            LOGGER.info("No lane with flowcell=%s, flowlane=%s, library=%s." % (flowcell, flowlane, lib.code))
            # Collect facility information
            try:
                facobj = Facility.objects.get(code=facility)
            except Facility.DoesNotExist, _err:
                raise SystemExit("No facility %s in repository." % facility)
            # Collect sequencing machine info
            if facility == 'EDG':
                try:
                    machine_obj = Machine.objects.get(code=instrument, platform='Illumina HiSeq X Ten')
                except Machine.DoesNotExist, _err:                    
                    LOGGER.info("Added Machine with code %s" % instrument)
                    machine_obj = Machine(code=instrument, platform='Illumina HiSeq X Ten', name='EdinburghX10')
                    machine_obj.save()
            else:
                try:
                    machine_obj = Machine.objects.get(code=instrument)
                except Machine.DoesNotExist, _err:                    
                    raise SystemExit("No machine %s in repository." % instrument)
            # Collect status object
            try:
                status_obj = Status.objects.get(code=status)
            except Status.DoesNotExist, _err:
                raise SystemExit("Status '%s' does not exist!\n" % status)
            # Ready to create lane
            lane = Lane(library  = lib,
                        flowcell = flowcell,                            
                        rundate  = rundate, # '2008-01-01', # Default. The value should
                        lanenum  = Lane.objects.next_lane_number(lib),
                        flowlane = int(flowlane),
                        paired = lib.paired,                        
                        genomicssampleid = sampleid, 
                        usersampleid = lib.code,
                        runnumber = runnumber,

                        facility = facobj,
                        seqsamplepf = '',
                        seqsamplebad = '',
                        
                        failed = False,
                        
                        machine = machine_obj,
                        status  = status_obj)
            lane.save()
            LOGGER.info("Created lane for flowcell=%s, flowlane=%s, library=%s." % (flowcell, flowlane, lib.code))
        return lane

    def date_to_sqldate(self, datestr):
        sqldate = ''
        if len(datestr) == 6:
            sqldate = '20%s-%s-%s' % (datestr[:2], datestr[2:4], datestr[4:6])
        return sqldate

    def ed_check(self):
        '''Downloads summary files for the project. Checks if all libraries have been registered in the repository,
        if files for any of the libraries have been downloaded or processed. Returns list of lines from all_md5sums.txt
        for which the files need to be downloaded and processed.'''
        
        # Download files in project root
        self._ed_parse_project_files()

        django.setup()
        
        # Check that libraries for all donumbers are in repository
        missing_lib=False
        for libcode in self.sample2library:
            # Find library in repository
            try:
                lib = Library.objects.search_by_name(self.sample2library[libcode])
            except Library.DoesNotExist, _err:
                missing_lib=True
                LOGGER.error("%s\tnoLibrary" % sample2library[libcode])
        if missing_lib:
            sys.exit(1)

        # For each file, check if lane exists (if not then create) and if data is already in repository or has been downloaded.
        for edstem in self.edfiles:            
            # self.ed_get_files_by_fastqfile(self.edfiles[edstem].file1, self.edfiles[edstem].file1_md5)

            # Parse one of the file names in pair and check if lane entry for it exists
            fname = self.edfiles[edstem].file1                
            (releasedate, sampleid, tmp, filename) = fname.split('/')
            (rundate, instrument, instrument_runnumber, flowcell_info, flowlane, adapterseq, tmp) = filename.split('_')
            parity = tmp.split('.')[0]
            if parity == 'R1':
                flowpair = '1'
            elif parity == 'R2':
                flowpair = '2'
            else:
                LOGGER.error('Failed to parse parity. Expected either \'R1\' or \'R2\', observed \'%s\'' % parity)
                sys.exit(1)                    
            rundate = self.date_to_sqldate(rundate)
            libcode = self.sample2library[sampleid]
            flowcell = flowcell_info[1:]
            flowcellpos = flowcell_info[:1]
            facility = 'EDG'
            status = 'new'
            notes = "ProjectName=%s; FlowcellPosition=%s; DataReleaseDate=%s" % (self.project, releasedate, flowcell[:1])

            # Find library in repository
            try:
                lib = Library.objects.search_by_name(libcode)
            except Library.DoesNotExist, _err:
                raise SystemExit("No library %s in repository." % libcode)
            # Find adapter by sequence
            try:
                adapter = Adapter.objects.get(sequence=adapterseq,protocol='TruSeq')
            except Adapter.DoesNotExist, _err:
                raise SystemExit("No TruSeq adapter with sequence \'%s\' in repository! Exiting." % adapterseq)
            if lib.adapter is None or lib.adapter.id != adapter.id:
                lib.adapter=adapter
                LOGGER.info("Adding adapter info to %s" % libcode)                    
            lib.paired=True
            lib.save()

            # Get (or register) lane
            lane = self.create_lane(lib, flowcell, flowlane, facility, status, rundate, sampleid, instrument_runnumber, instrument, notes)
            self.edfiles[edstem].laneid = lane.id
            self.edfiles[edstem].status = lane.status.code

            # Create repository file stem
            repfn = "%s_%s%02d" % (lib.filename_tag,
                                        facility,
                                        lane.lanenum)
            #Just replace spaces for now as e.g. UCSC upload fails in
            # these cases. Also forward slashes. And parentheses/semicolons.
            sanity_re = re.compile(r'([ \/\(\);]+)')
            repstem     = sanity_re.sub('_', repfn)
            # Save repository filename for further processing        
            self.edfiles[edstem].repstem = repstem

            LOGGER.info("%s\t%s\t%d\t%s", edstem, lane.status.code, lane.id, repstem)

    def ed_download(self, print_download_commands_only=False):
        '''Downloads all project related files.'''

        # Set some variables for submitting the download jobs to cluster.
        submitter = ClusterJobSubmitter()

        # Install information about files locations, ids etc.
        self.ed_check()
        
        # Set some variables for threaded downloading.
        # We will try to control for the number of threads by
        # setting new download jobs to depend on complete of previous download jobs
        # The assumption is that in average the run time of all download threads will be the same.
        tnr = 0
        jobids = []
        newids = []
            
        for edstem in self.edfiles:
            
            status = self.edfiles[edstem].status
                    
            if status != 'new':
                if status == 'downloaded':
                    LOGGER.info("Skipping %s (%s)", edstem, status)
                if status == 'processed':
                    LOGGER.info("Skipping %s (%s)", edstem, status)
                continue

            ## Prepare download command for for submision to the cluster
            cmd = 'cs_edinburgh_download.py -a --file1 %s --file2 %s --file1_md5 %s --file2_md5 %s -p %s -l %d' % (self.edfiles[edstem].file1, self.edfiles[edstem].file2, self.edfiles[edstem].file1_md5, self.edfiles[edstem].file2_md5, self.project, self.edfiles[edstem].laneid)
            if print_download_commands_only:
                print command
                continue
            # Submit download job
            if jobids:
                jobid = submitter.submit_command(cmd=cmd, mem=1000, auto_requeue=False, depend_jobs=[jobids[tnr]])
            else:
                jobid = submitter.submit_command(cmd=cmd, mem=1000, auto_requeue=False)
            LOGGER.info("Setting up downloads for %s ... (jobid=%s)", edstem, jobid)
            newids.append(int(jobid))
            tnr += 1
            if tnr == self.athreads:
                jobids = newids
                newids = []
                tnr = 0

    def ed_process(self, print_commands_only=False):        
        '''Process downloaded files'''
        # NB! In future, the commands in this function sould be added to the end of function 'ed_get_files_by_fastqfile'
        
        # Set some variables for submitting the download jobs to cluster.
        submitter = ClusterJobSubmitter()

        # Check the project files first
        self.ed_check()

        django.setup()

        # Fetch archive location in the file system.
        try:
            archive = ArchiveLocation.objects.get(name=self.aname)
        except ArchiveLocation.DoesNotExist, _err:
            raise SystemExit("Archive location \'%s\' does not exist!" % self.aname)
            
        # Set some variables for threaded downloading.
        # We will try to control for the number of threads by
        # setting new download jobs to depend on complete of previous download jobs
        # The assumption is that in average the run time of all download threads will be the same.
        tnr = 0
        jobids = []
        newids = []

        for edstem in self.edfiles:

            # Process only data that has been downloaded
            if self.edfiles[edstem].status != 'downloaded':
                continue

            # Build full path to repository file names to be created
            repstem = self.edfiles[edstem].repstem
            rfastq1 = os.path.join(self.destination, repstem + 'p1.fq')
            rfastq2 = os.path.join(self.destination, repstem + 'p2.fq')
            rfastq1_md5 = os.path.join(self.destination, rfastq1 + '.md5')
            rfastq2_md5 = os.path.join(self.destination, rfastq2 + '.md5')
            summaryFile = os.path.join(self.destination, rfastq1 + '.summary')

            # Build full path to downloaded files
            (releasedate, sampleid, tmp, filename1) = self.edfiles[edstem].file1.split('/')
            (releasedate, sampleid, tmp, filename2) = self.edfiles[edstem].file2.split('/')                        
            dfile1 = os.path.join(self.destination, filename1)
            dfile2 = os.path.join(self.destination, filename2)
        
            # Set up directory for archiving the file
            libcode = repstem.split('_')[0]
            apath = os.path.join(archive.root_path, libcode)
            if not os.path.exists(apath):
                os.makedirs(apath)
                LOGGER.info("Creating archive directory for %s" % apath)
            # Build full path archived files
            afastq1 = os.path.join(apath, os.path.split(rfastq1)[1] + '.gz')
            afastq2 = os.path.join(apath, os.path.split(rfastq2)[1] + '.gz')

            # Build commands
            # Uncompress and compute md5 sum on fly. Move compressed file to archive. Delete .md5 and .done files generated during download.
            cmd = 'time zcat %s | tee %s | md5sum > %s && time mv %s %s && rm %s' % (dfile1, rfastq1, rfastq1_md5, dfile1, afastq1, dfile1 + '.*')
            cmd += ' && time zcat %s | tee %s | md5sum > %s && time mv %s %s && rm %s' % (dfile2, rfastq2, rfastq2_md5, dfile2, afastq2, dfile2 + '.*')

            # Compute summary for file2
            cmd += ' && time cat %s | summarizeFile > %s' % (rfastq1, summaryFile)

            # Register processed files in repository
            cmd += ' && time cs_addFile.py -m -s %s --archive ebiark %s %s && rm %s' % (summaryFile, rfastq1, rfastq2, summaryFile)

            # Set lane status complete
            cmd += ' && time communicateStatus.py --laneid %d --status complete' % self.edfiles[edstem].laneid
        
            if print_commands_only:
                print "%s" % cmd
                continue

            # Submit command to cluster
            if jobids:
                jobid = submitter.submit_command(cmd=cmd, mem=1000, auto_requeue=False, depend_jobs=[jobids[tnr]], mincpus=1)
            else:
                jobid = submitter.submit_command(cmd=cmd, mem=1000, auto_requeue=False, mincpus=1)
            LOGGER.info("Submitting job with commands \'%s\' (jobid=%s)" % (cmd, jobid))
            newids.append(int(jobid))
            tnr += 1
            if tnr == self.athreads:
                jobids = newids
                newids = []
                tnr = 0

    def ed_get_files_by_fastqfile(self, file1, file1_md5, laneid, file2=None, file2_md5=None):
        '''Downloads all files for the userid'''
        # For each fastq file in Edinburgh Genomics, there is also a *_R1_fastqc.html file to download

        file1_fastqc = file1.rstrip('.fastq.gz') + '_fastqc.html'
        rfiles = [file1, file1_fastqc]
        rmd5s = [file1_md5, None]

        fail_command = 'cs_edinburgh_download.py -a --file1 %s --file1_md5 %s -p %s -l %d' % (file1, file1_md5, self.project, laneid)

        if file2 is not None:
            file2_fastqc = file2.rstrip('.fastq.gz') + '_fastqc.html'
            rfiles = rfiles + [file2, file2_fastqc]
            rmd5s = rmd5s + [file2_md5, None]

            fail_command += " --file2 %s --file2_md5 %s" % (file2, file2_md5)

        failed = False        

        # Download files in rfiles
        for rfpath, md5sum in zip(rfiles, rmd5s):
            attempts = 0
            while attempts < self.maxattempts:
                # Depending if the file to be downloaded is .html skip or download also associated .md5 file
                if rfpath.endswith('.html'):                    
                    res = self.ed_get_file_with_md5(rfpath, md5sum, download_md5=False)
                else:
                    res = self.ed_get_file_with_md5(rfpath, md5sum, download_md5=True)
                if res == 0:
                    break
                else:            
                    attempts += 1
                    if attempts == self.maxattempts:
                        LOGGER.error("Giving up (after %d attempts) on trying to download %s" % rpath)

                        # Set download failed
                        cmd = 'communicateStatus.py --laneid %d --status downloading_failed' % laneid
                        # Record failed command (in case its of any use)
                        cmd += ' && echo %s >> %s' % (fail_command, self.failedcommands)
                        run_command(cmd, shell=True)
                        failed = True
                    else:
                        LOGGER.error("Trying to download %s again." % rpath)
                        
        # Register files as downloaded
        if not failed:
            # Get flowcell and machine info from read header of fastq file
            # (machine, flowcell, flowlane) = parse_read_header(fqfile1)
            # 5. Set flowcell status 'downloaded' or 'failed download'
            # cmd = 'communicateStatus.py --status complete --flowcell %s --flowlane 0 --library %s --facility EDG' % (flowcell, userid)
            cmd = 'communicateStatus.py --laneid %d --status downloaded' % laneid
            call_subprocess(cmd, shell=True)
            
    def ed_get_file_with_md5(self, rfile, md5sum, download_md5=True, validate_md5=True):
        '''Initiates download for remote file rfile, rfname.md5 for the project. Returns 0 in case downloaded fname matches downloaded md5.'''

        # NB! Note that the code below is agnostic to whether dowloading of the file or md5 file may have failed.        
        ret = 0
        
        # build full file path for remote and local files
        rfpath = os.path.join(self.project, rfile)
        fname = os.path.split('/' + rfile)[1]
        lfpath = os.path.join(self.destination, fname)

        if os.path.exists(lfpath + '.done') and os.path.exists(lfpath):
            LOGGER.info("Skipping download for %s. File already exists." % fname)
        else:
            # if either of the files exists, remove before downloading again.
            if os.path.exists(lfpath + '.done'):
                os.unlink(lfpath + '.done')
            if os.path.exists(lfpath):
                os.unlink(lfpath)
            ret = self.ed_download_file(rfpath, lfpath)
            if ret > 0:
                return ret

        if download_md5:
            # download md5 for the file
            fnmd5 = lfpath + '.md5'
            if not os.path.exists(fnmd5):
                ret = self.ed_download_file(rfpath + '.md5', fnmd5)
                if ret > 0:
                    return ret
            # assuming that if md5 has been downloaded in past, the fastq has also been checked.
            else:
                LOGGER.info("Skipping download for %s. File already exists. Assuming md5 sums have been compared." % fname)
                return 0
            
            if validate_md5:
                # compute md5 for the downloaded file
                md5 = compute_md5(lfpath)                
                with open(fnmd5, 'rb') as fh:
                    for line in fh:
                        line = line.rstrip('\n')
                        (md5col, fcol) = line.split()
                        if md5col == md5:
                            if md5sum is not None and md5sum == md5:                            
                                LOGGER.info("MD5 is correct!")
                                ret = 0
                            else:
                                LOGGER.error("Md5 sum computed for file (%s) and provided on command line (%s) are different!", md5, md5sum)
                                ret = 1
                        else:
                            LOGGER.error("Md5 sum computed for file (%s) and provided on command line (%s) are different!", md5, md5sum)
                            ret = 1
                            break
                if ret == 1:
                    LOGGER.error("Removing corrupt file %s", lfpath)
                    os.unlink(lfpath)
                    os.unlink(lfpath + '.done')
                    os.unlink(fnmd5)
                    os.unlink(fnmd5 + '.done')
        return ret

    def ed_download_file(self, rfpath, lfpath):
        '''Downloads file with name fname'''
        
        LOGGER.info("Downloading file %s" % rfpath)
        acmd = "ascp -q -T -p -P %s -O %s -l %s %s@%s:%s %s" % (self.aPport, self.aOport, self.arate, self.credentials['username'], self.ahost, rfpath, lfpath)
        (stdout, stderr, retcode) = run_command(acmd, shell=True, path=self.path)
        # call_subprocess(acmd, shell=True, path=self.path) # Make sure ASPERA_SCP_PASS is available in shell.     
        if retcode == 0:
            LOGGER.info("Download complete. File location: %s" % lfpath )
            # Mark file on disk as downloaded successfully
            cmd = 'touch %s.done' % lfpath
            LOGGER.info("Touching file \'%s\'" % cmd )
            (stdout, stderr, retcode) = run_command(cmd, shell=True, path=self.path)
        else:
            LOGGER.error("Download failed! Could it be that file %s does not exist in the server?" % rfpath)

        return retcode


##################  M A I N   P R O G R A M  ######################

if __name__ == '__main__':

    from argparse import ArgumentParser
    
    PARSER = ArgumentParser(description='Submission of project related files to ENA')
    PARSER.add_argument('-p', '--project', dest='project', type=str, help='Name of the project (subdirectory) in Edinburgh aspera server where the library files are located.', required=True, default='')
    PARSER.add_argument('--file1', dest='file1', type=str, help='Download file and file1.md5 in the project subdirectory. A file.done is created in file system on success. NB! This variable expects a full path to the file within the project subdirectory!')
    PARSER.add_argument('--file2', dest='file2', type=str, help='Download file and file2.md5 in the project subdirectory. A file.done is created in file system on success. NB! This variable expects a full path to the file within the project subdirectory!')
    PARSER.add_argument('-a', '--all-for-file', dest='allforfile', action='store_true', help='Downloads all files for a fastq file (fastq itself, its md5 and fastqc report.', default=False)
    PARSER.add_argument('-l', '--lane_id', dest='laneid', type=int, help='Lane id associated with the file.', default=None)
    PARSER.add_argument('-s', '--single_file', dest='single_file', action='store_true', help='Download file in the project subdirectory.', default=False)
    PARSER.add_argument('--file1_md5', dest='file1_md5', type=str, help='Md5 sum for file1.', default=None)
    PARSER.add_argument('--file2_md5', dest='file2_md5', type=str, help='Md5 sum for file1.', default=None)
    PARSER.add_argument('-A', '--archive-name', dest='aname', type=str, help='Arhcive name where to save the fastq files.', default='ebiark')
    PARSER.add_argument('--test', dest='print_commands_only', action='store_true', help='Only print the download commands of fastq files.', default=False)
    PARSER.add_argument('--process', dest='process', action='store_true', help='Process downloaded files.', default=False)

    PARSER.add_argument('--check', dest='check', action='store_true', help='Check if libraries for all files part of the project exist in repository. Where necessary, create new lanes and set the status \'new\'.', default=False)
    PARSER.add_argument('--download', dest='download', action='store_true', help='Download all files for which the lanes are marked \'new\'.', default=False)
    PARSER.add_argument('-t', '--threads', dest='threads', type=int, help='Number of threads. Overrides default (athreads) from configuration.', default=10)

    ARGS = PARSER.parse_args()
    edd = ed_data_handler(ARGS.project, threads=ARGS.threads, aname=ARGS.aname)

    # Check and register files that need to be downloaded
    if ARGS.check:
        edd.ed_check()
        sys.exit(0)

    # Download files to incoming
    if ARGS.download:
        edd.ed_download(print_download_commands_only=ARGS.print_commands_only)
        sys.exit(0)

    # Process files successfully downloaded to incoming.
    if ARGS.process:
        edd.ed_process(print_commands_only=ARGS.print_commands_only)
        sys.exit(0)

    # Download for only one fastq file
    if ARGS.file1:
        if ARGS.allforfile:
            if ARGS.laneid is None:
                sys.stderr.write("Lane id for for the file missing! Use -l or --lane_id option! Exiting!")
                sys.exit(1)
            else:
                edd.ed_get_files_by_fastqfile(ARGS.file1, ARGS.file1_md5, ARGS.laneid, file2=ARGS.file2, file2_md5=ARGS.file2_md5)
        elif ARGS.single_file:            
            edd.ed_get_file_with_md5(ARGS.file1, ARGS.file1_md5, download_md5=False)
        else:
            edd.ed_get_file_with_md5(ARGS.file1, ARGS.file1_md5)


