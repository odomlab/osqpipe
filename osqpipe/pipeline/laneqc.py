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

'''Classes and functions used to generate QC reports for data
associated with sequencing lanes.'''

import os
import re
import tarfile

from tempfile import mkdtemp
from shutil import rmtree, move, copy
from pkg_resources import Requirement, resource_filename

from django.db import transaction, models
from ..models import Program, LaneQC, QCfile, Filetype, Lanefile, DataProvenance, Datafile, DataProcess
from osqutil.progsum import ProgramSummary
from osqutil.utilities import checksum_file, call_subprocess, rezip_file, set_file_permissions, transfer_file
from osqutil.config import Config
from osqutil.setup_logs import configure_logging

CONFIG = Config()
LOGGER = configure_logging('laneqc')

class QCReport(object):
  '''
  Abstract superclass handling all QC reports. Note that this is
  implemented as a context manager, so you would typically use
  subclasses in the following way::

    with LaneFastQCReport(target=l, program_name='fastqc') as rep:
      rep.insert_into_repository()
  '''

  __slots__ = ('target', 'workdir', 'output_files', 'program_name', 'path',
               'program_params', '_dbprog', '_delete_workdir','output_md5s','move_files')

  data_process     = DataProcess # for the benefit of pylint
  target_name      = None
  data_file        = Datafile    # for the benefit of pylint
  file_target_name = None
  
  def __init__(self, target, program_name, path=None, program_params='',
               workdir=None, move_files=True):

    self.target         = target
    self.program_name   = program_name
    self.program_params = program_params
    self.path           = path
    
    self.output_files   = []
    self.output_md5s    = []

    self.move_files = move_files
    
    self.workdir = workdir
    if workdir is not None:
      self._delete_workdir = False
    else:
      if self.move_files == False:
        raise StandardError("Not moving files from temporary directory to be deleted does not make sense!")        

    # This checks that the specified program exists, and where it
    # yields some kind of meaningful version info will record that.
    progdata = ProgramSummary(program_name, path=path)
    
    # This is a little vulnerable to correct version parsing by
    # progsum.
    try:
      self._dbprog = Program.objects.get(program = progdata.program,
                                         version = progdata.version,
                                         current = True)
    except Program.DoesNotExist, _err:
      raise StandardError(("Unable to find current %s program (version %s)"
                           + " record in the repository")
                          % (progdata.program, progdata.version))
      
  def __enter__(self):
    if self.workdir is None:
      self.workdir = mkdtemp(dir=CONFIG.tmpdir)
      LOGGER.debug("Working directory is %s", self.workdir)
      self._delete_workdir = True
    return self

  def generate(self):
    '''
    Generate report and insert output files into self.output_files. To
    be implemented in subclasses.
    '''
    raise NotImplementedError()

  @transaction.atomic
  def insert_into_repository(self, move_files=True):
    '''Insert self.output_files into the database.'''

    if len(self.output_files) == 0:
      self.generate()

    params = { self.target_name : self.target }
    qcobj  = self.data_process.objects.create(**params)
    DataProvenance.objects.create(program      = self._dbprog,
                                  parameters   = self.program_params,
                                  rank_index   = 1,
                                  data_process = qcobj)

    for i in range(len(self.output_files)):
      fname = self.output_files[i]
      if len(self.output_md5s) != len(self.output_files):
        checksum = None
      else:
        checksum = self.output_md5s[i]
        
      LOGGER.info("Inserting %s", fname)
      # Note: this will fail if multiple types match.
      ftype = Filetype.objects.guess_type(fname)

      if os.path.isabs(fname):
        fpath = fname
      else:
        fpath = os.path.join( self.workdir, fname )
        
      if checksum is None or checksum == '':
        checksum = checksum_file(fpath)

      fparms = { self.file_target_name : qcobj,
                 'filename'            : os.path.split(fname)[1],
                 'checksum'            : checksum,
                 'filetype'            : ftype }
      fobj = self.data_file(**fparms)

      fobj.save()

      if move_files:      
        # Zip up the file if necessary.
        if ftype.gzip and os.path.splitext(fname)[1] != CONFIG.gzsuffix:
          fpath = rezip_file(fpath)
        if self.move_files:
          dest    = fobj.repository_file_path
          # destdir = os.path.dirname(dest)
          # if not os.path.exists(destdir):
          #    os.makedirs(destdir)
          # move(fpath, dest)
          # set_file_permissions(CONFIG.group, dest)
          if os.path.isabs(dest):
            dest = os.path.split(dest)[0] + '/'
          transfer_file(fpath, "%s@%s:%s" % (CONFIG.user, CONFIG.datahost, dest), set_ownership=True) # note that transfer_file sets destination file permissions as in CONF

  def __exit__(self, exctype, excvalue, traceback):
    if self._delete_workdir:
      try:
        LOGGER.debug("Deleting working directory %s", self.workdir)
        rmtree(self.workdir)
      except Exception, _err:
        LOGGER.warning("Unable to delete working directory %s", self.workdir)

class LaneQCReport(QCReport):
  '''
  Abstract class handling all lane-based QC Reports.
  '''
  data_process     = LaneQC
  target_name      = 'lane'
  data_file        = QCfile
  file_target_name = 'laneqc'

class LaneFastQCReport(LaneQCReport):
  '''
  Concrete LaneQCReport subclass implementing fastqc report
  generation. See the superclass for usage notes.
  '''
  def __init__(self, fastqs=None, program_name='fastqc', *args, **kwargs):
    '''
    The fastqs attribute is to allow callers to override Lane objects
    which have no fastq files attached; for example, when loading
    external data from public repositories such as GEO or
    ArrayExpress, we want to store the FastQC report but not the
    originial fastq files. In such cases the fastq files should be
    passed in alongside the Lane object.
    '''
    super(LaneFastQCReport, self).\
        __init__(program_name=program_name, *args, **kwargs)
    self.fastqs = fastqs

  def generate(self):

    if self.fastqs is not None:
      LOGGER.debug('Using the provided fastq files.')
      fns = self.fastqs
    else:
      LOGGER.debug('Using the fastq files stored in the repository.')
      lanefiles = Lanefile.objects.filter(lane=self.target,
                                          filetype__code='fq')
      assert(len(lanefiles) > 0)
      fns = [ x.repository_file_path for x in lanefiles ]

    self.run_fastqc(fns)
    self.postprocess_results(fns)

  def run_fastqc(self, fns, threads=2):

    """Executes fastqc report generation."""

    assert(len(fns) > 0)

    if len(fns) < threads:
      threads = len(fns)
    cmd = [ self.program_name,
            '-q',
            '-t', threads,
            '-o', self.workdir ]

    if len(self.program_params) > 0:
      cmd.extend( self.program_params.split() )

    cmd.extend(fns)

    cmd = [ str(x) for x in cmd ]
    LOGGER.info("Running FastQC command: %s", " ".join(cmd))
    call_subprocess(cmd, path=self.path)

  def postprocess_results(self, fns):

    '''Checks for output, add to self.output_files. Note that we want
    the compressed archive to be gzipped (*.tar.gz), not zipped. We
    also want the fastqc_report.txt file stored separately and
    uncompressed.'''

    for fpath in fns:

      fname = os.path.split(fpath)[1]
      fname = re.sub(r'\.gz$', '', fname)

      # FastQC strips '.fastq' but not '.fq', so we only remove the former here.
      fname = re.sub(r'\.fastq$', '', fname)

      base  =  "%s_fastqc" % fname
      bpath = os.path.join(self.workdir, base)

      if not os.path.exists(bpath):
        raise StandardError("Expected output directory not found: %s" % bpath)

      # Sort out the tar-gzipped archive.
      gzarch = "%s.tar" % bpath
      tar    = tarfile.open(gzarch, mode='w')

      # A little jimmying around so we only get the directory we want.
      pwd = os.getcwd()
      os.chdir(self.workdir)
      tar.add(base)
      os.chdir(pwd)

      tar.close()
      self.output_files.append(gzarch)
      self.output_md5s.append(checksum_file(gzarch))

      # The text file containing summary results. Useful for analyses.
      resfile = "%s.txt" % bpath
      copy(os.path.join(bpath, 'fastqc_data.txt'), resfile)
      self.output_files.append(resfile)
      self.output_md5s.append(checksum_file(resfile))
      
      # Generating a PDF for our end-users.
      html = os.path.join(bpath, 'fastqc_report.html')
      pdf  = "%s.pdf" % bpath

      # FIXME resource_filename is a little brittle, would
      # resource_string be better?
      cmd = [ 'wkhtmltopdf-amd64',
              '--user-style-sheet',
              resource_filename(Requirement.parse('osqpipe'),
                                'osqpipe/pipeline/fastqc_pdf_styles.css'),
              html, pdf ]
      call_subprocess(cmd, path=self.path)
      self.output_files.append(pdf)
      self.output_md5s.append(checksum_file(pdf))
