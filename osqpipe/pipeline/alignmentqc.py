'''Classes and functions used to generate QC reports for data
associated with sequencing lanes.'''

import os
import re

from tempfile import NamedTemporaryFile
from shutil import move

from django.db import transaction
from ..models import Program, AlignmentQC, AlnQCfile, Filetype, Alnfile, DataProvenance
from .laneqc import QCReport
from osqutil.progsum import ProgramSummary
from osqutil.utilities import checksum_file, call_subprocess, rezip_file, set_file_permissions
from osqutil.config import Config
from osqutil.setup_logs import configure_logging

CONFIG = Config()
LOGGER = configure_logging('alignmentqc')

class AlignmentQCReport(QCReport):
  '''
  Abstract class handling all alignment-based QC Reports.
  '''
  data_process     = AlignmentQC
  target_name      = 'alignment'
  data_file        = AlnQCfile
  file_target_name = 'alignmentqc'

class AlignmentCrossCorrReport(AlignmentQCReport):
  '''
  Concrete AlignmentQCReport subclass implementing cross-correlation
  analysis plots. See the superclass for usage notes.
  '''
  def __init__(self, bams=None, program_name='run_spp_nodups.sh', *args, **kwargs):
    '''
    The bams attribute is to allow callers to override Alignment
    objects which have no fastq files attached; for example, when
    loading external data from public repositories such as GEO or
    ArrayExpress, we want to store the cross-correlation report but
    not the originial bam file. In such cases the bam files should be
    passed in alongside the Alignment object.
    '''
    super(AlignmentCrossCorrReport, self).\
        __init__(program_name=program_name, *args, **kwargs)
    self.bams = bams

  def generate(self):

    if self.bams is not None:
      LOGGER.debug('Using the provided bam file.')
      fns = self.bams
    else:
      LOGGER.debug('Using the bam files stored in the repository.')
      alnfiles = Alnfile.objects.filter(alignment=self.target,
                                        filetype__code='bam')
      assert(len(alnfiles) > 0)
      fns = [ x.repository_file_path for x in alnfiles ]

    self.run_analysis(fns)

  def run_analysis(self, fns):

    """Executes cross-correlation report generation."""

    assert(len(fns) == 1)

    basefn = os.path.splitext(fns[0])[0]
    pdf    = "%s_xcor.pdf" % basefn
    out    = "%s_xcor.txt" % basefn

    with NamedTemporaryFile(suffix='.bam', dir=self.workdir) as tempbam:
      with NamedTemporaryFile(suffix='.txt', dir=self.workdir) as tempout:

        mdcmd = [ 'picard',
                  'MarkDuplicates',
                  'I=%s' % fns[0],
                  'O=%s' % tempbam.name,
                  'M=%s' % tempout.name,
                  'REMOVE_DUPLICATES=true',
                  'VALIDATION_STRINGENCY=SILENT' ]

        LOGGER.debug('Removing bam file duplicate reads: %s', fns[0])

        call_subprocess(mdcmd, path=self.path)

      # tempout closes and is deleted.
      cmd = [ self.program_name,
              '-c=%s' % tempbam.name,
              '-savp=%s' % pdf,
              '-out=%s' % out ]

      if len(self.program_params) > 0:
        cmd.extend( self.program_params.split() )

      cmd = [ str(x) for x in cmd ]
      LOGGER.info("Running Cross-correlation analysis: %s", " ".join(cmd))
      call_subprocess(cmd, path=self.path)

    # tempbam closes and is deleted.
    self.output_files.extend([out, pdf])
