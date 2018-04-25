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

'''
A script designed to automate the process of creating axtNet files for
a LASTZ-based pairwise alignment between two genomes. See this page
for the basic model:

http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto

This script is heavily reliant on the kent src utils for working with
fasta and 2bit files, and obviously the lastz binary itself.
'''

import os
import sys
import socket
import re

from time import sleep
from hashlib import md5

# We can hook into our main pipeline code for our external binary and
# cluster process handling. This also takes care of cluster host
# configuration for us.

# Required binaries in the cluster PATH: lastz, scp, ssh, rm, bsub, touch, find
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from osqutil.utilities import call_subprocess, bash_quote, \
  CalledProcessError
from osqutil.setup_logs import configure_logging
from osqpipe.pipeline.bwa_runner import ClusterJobManager

from tempfile import gettempdir
from logging import INFO, DEBUG

LOGGER = configure_logging(level=INFO)

MASKTAG = '_trfBig_masked'

def filebasename(fname):
  '''
  Returns the base name for a file with its extension removed.
  '''
  return os.path.splitext(os.path.basename(fname))[0]

def splitpath(path):
  '''
  Recursive function to do what I personally think os.path.split()
  should do.
  '''
  parts = os.path.split(path)
  if parts[0] in ('.', '/', ''):
    return parts[1:]
  else:
    return splitpath(parts[0]) + parts[1:]

class LastzAligner(ClusterJobManager):
  '''
  Class to handle all the steps required for generating an axt-format
  net alignment file using lastz as the aligner.
  '''
  # local_tempdir will need to be able to handle around 75GB when aligning
  # two typical mammalian genomes.

  def __init__(self, from_genome, to_genome, hsp_thresh=3000,
               length_limit=None, linear_gap='loose', local_tempdir=None,
               resume=False, *args, **kwargs):

    super(LastzAligner, self).__init__(*args, **kwargs)

    self.from_genome  = from_genome
    self.to_genome    = to_genome
    self.hsp_thresh   = hsp_thresh
    self.length_limit = length_limit
    self.linear_gap   = linear_gap

    # Flag used to tell the object to fill in missing lav files by
    # resubmitting to the cluster, rather than just working with
    # what's available.
    self.resume = resume

    systempdir = gettempdir() if local_tempdir is None else local_tempdir
    self.local_tempdir = os.path.join(systempdir, str(os.getpid()))
    os.mkdir(self.local_tempdir) # Fails on pre-existing directory.

    self._tempfiles  = []
    self._chr_sizes  = {}

  def split_chrs(self, fasta, dryrun=False):
    '''
    Split a designated fasta file by chromosome. Returns a list of the
    generated fasta files. Any chromosome whose sequence exceeds
    self.length_limit will be split appropriately. Calling with
    dryrun=True returns a list of files which would have been created;
    this may be useful when deciding on an appropriate length_limit
    parameter.
    '''
    LOGGER.info("Splitting fasta by chromosome: %s", fasta)
    # N.B. the trailing '/' is important here:
    wdir = os.path.join(self.local_tempdir,
                        '%s_chr_split/' % filebasename(fasta))
    if not dryrun:
      os.mkdir(wdir) # Fails on pre-existing directory.
    self._tempfiles.append(wdir)

    outfiles = []
    handle = open(fasta, 'rU')
    for chromosome in SeqIO.parse(handle, 'fasta'):

      # Check whether we need to split the chromosome.
      seqlen = len(chromosome.seq)
      if self.length_limit and seqlen > self.length_limit:

        # Figure out how many chunks we need.
        denom = 2
        while (float(seqlen)/denom) > self.length_limit:
          denom += 1

        # Output the sequences
        for segnum in range(denom):
          start = (segnum * (seqlen/denom)) + 1
          end   = min(seqlen, (segnum+1) * (seqlen/denom))

          # This filename format will be parsed later, in
          # process_lavs_to_psl. The filename coordinate needs to be
          # added to the output psl coords.
          new_id  = "%s_+%d" % (chromosome.id, start-1)
          chrfile = os.path.join(wdir, "%s.fa" % new_id)
          chrseg  = chromosome[start-1:end]
          chrseg.id = new_id
          if not dryrun:
            with open(chrfile, 'w') as chrfh:
              SeqIO.write([chrseg], chrfh, 'fasta')
          outfiles.append(chrfile)

      else:

        # If chromosome is small enough, just dump it out in a single file.
        chrfile = os.path.join(wdir, "%s.fa" % chromosome.id)
        if not dryrun:
          with open(chrfile, 'w') as chrfh:
            SeqIO.write([chromosome], chrfh, 'fasta')
        outfiles.append(chrfile)
        
    return outfiles

  def mask_tandem_repeats(self, fasta):
    '''
    Runs trfBig over the designated fasta file. Should return the
    newly-generated masked fasta file name. Runs quite slowly, so we
    keep the outputs following 2bit conversion.
    '''
    LOGGER.info("Masking tandem repeats for fasta: %s", fasta)
    curdir = os.getcwd()

    # trfBig writes to current working directory a lot.
    os.chdir(self.local_tempdir)
    maskfn = os.path.splitext(fasta)[0] + MASKTAG + '.fa'
    cmd    = [ 'trfBig', fasta, maskfn ]
    call_subprocess(cmd, tmpdir=self.local_tempdir, path=os.environ['PATH'])
    os.chdir(curdir)
    return maskfn

  def convert_to_2bit(self, fasta, workdir=None):
    '''
    Runs faToTwoBit on the designated fasta file; returns the name of
    the output 2bit file.
    '''
    LOGGER.info("Converting fasta to 2bit: %s", fasta)
    if workdir is None:
      workdir = self.local_workdir

    twobitfn = os.path.join(workdir, filebasename(fasta) + '.2bit')
    cmd  = [ 'faToTwoBit', '-noMask', fasta, twobitfn ]
    call_subprocess(cmd, tmpdir=self.local_tempdir, path=os.environ['PATH'])
    return twobitfn

  def make_cluster_filename(self, localfile):
    '''
    Generate a unique filename to be used on the cluster, without
    unnecessarily divulging local file paths.
    '''
    pathbits = os.path.split(localfile)
    hasher   = md5()
    hasher.update(pathbits[0])
    clusterfile = "%d_%s_%s" % (os.getpid(), hasher.hexdigest(), pathbits[1])
    return clusterfile

  def align(self, from_list, to_list, omit_list=None):
    '''
    Actually run the alignment. Requires lastz on the cluster, and of
    course bsub/LSF et al. Note that we will generate one lastz
    process for each chr-chr combination, so this should lend itself
    well to a clustered solution.
    '''
    # Just in case, for convenience.
    if type(from_list) in (str, unicode):
      from_list = [from_list]
    if type(to_list) in (str, unicode):
      to_list = [to_list]
    if omit_list is None:
      omit_list = []

    # Make sure the filenames on the cluster won't easily collide.
    cluster_from = [ self.make_cluster_filename(x) for x in from_list ]
    cluster_to   = [ self.make_cluster_filename(x) for x in to_list   ]

    # Note that we copy all the files to the cluster even if we only
    # want to repeat a handful of alignments; managing the files is
    # simpler that way.
    LOGGER.info("Copying files to cluster server.")
    self.submitter.remote_copy_files(filenames=from_list + to_list,
                                   destnames=cluster_from + cluster_to)

    job_ids  = []
    lavfiles = []
    for from_num in range(len(from_list)):
      for to_num in range(len(to_list)):

        # Files on localhost
        from_file = from_list[from_num]
        to_file   = to_list[to_num]

        # Files on the cluster
        from_clust = cluster_from[from_num]
        to_clust   = cluster_to[to_num]

        outfile = "%s_%s.lav" % (filebasename(from_file),
                                 filebasename(to_file))

        if outfile in omit_list:
          LOGGER.warning("Skipping pre-existing lav file %s...", outfile)
          lavfiles.append(outfile)
          continue

        ## FIXME consider the --inner option here (ensembl-compara
        ## appears to use --inner=2200).
        LOGGER.info("Launching alignment (%s : %s).",
                    from_file, to_file)
        clusterout = "%d_%s" % (os.getpid(), outfile)

        ## We use this file to monitor lastz completion, to
        ## disambiguate lastz failure from scp failure. FIXME if this
        ## turns out to be scp failure we can add a final re-try to
        ## the monitor job.
        clusterdone = clusterout + '.done'

        ## Note that using --chain here appears to be undesirable
        ## since the lastz chaining implementation is rather too
        ## simplistic for our purposes (see lastz docs).
        cmd = ['lastz',
               to_clust, from_clust, # This is the correct order.
               '--format=lav',
               '--hspthresh=%d' % self.hsp_thresh,
               '--output=%s'    % clusterout]

        sshcmd = self.return_file_to_localhost(clusterout, outfile,
                                               execute=False)
        LOGGER.debug(sshcmd)
        cmd = " ".join(cmd) + (' && touch %s && %s && rm %s %s'
                               % (clusterdone, sshcmd, clusterout, clusterdone))

        # 4GB is the default max mem for lastz. Setting mem=4000 means
        # some larger alignments fail silently; using 5000 seems much
        # more robust on our cluster.
        job_ids.append(self.submitter.submit_command(cmd=cmd,
                                                     mem=self.memsize * 1024,
                                                     auto_requeue=False,
                                                     time_limit=self.time_limit))
        lavfiles.append(outfile)

        # Reduce the rate of cluster job submission, if desired.
        sleep(self.throttle)

    # Caller code tends to assume these paths are absolute.
    lavfiles = [ os.path.join(self.local_workdir, x) for x in lavfiles ]
        
    return (job_ids, lavfiles, cluster_from + cluster_to)

  def convert_to_psl(self, lav):
    '''
    Converts an input lav file to a temporary psl file.
    '''
    pslfn = os.path.join(self.local_tempdir, filebasename(lav) + '.psl')
    cmd = [ 'lavToPsl', lav, pslfn ]
    call_subprocess(cmd, tmpdir=self.local_tempdir, path=os.environ['PATH'])
    return pslfn # Delete this file in the caller code.

  def process_lavs_to_psl(self, lavs):
    '''
    Convert .lav files to .psl, swap query and target, and split on
    target.
    '''
    # Convert lav files to psl, concatenate.
    LOGGER.info("Reorganising lav files into psl files.")
    psls = [ self.convert_to_psl(x) for x in lavs ]
    allpsl = os.path.join(self.local_tempdir, 'all.psl')

    # Concatenate the files. We take this opportunity to strip out the
    # junk we've added to the chromosome names.
    def repl(match):
      '''
      Regex replace function.
      '''
      return "\t%s\t" % match.group(1)

    from_sizes = self.get_chr_sizes_dict(self.from_genome)
    to_sizes   = self.get_chr_sizes_dict(self.to_genome)

    # We allow for any pid prefix so we can restart in a new process
    # if needed. Also allow for genome/chrN_trfBig_masked to support
    # fill-in files generated locally.
    genstr   = (r'(?:%s|%s)' % (filebasename(self.from_genome),
                                filebasename(self.to_genome)))
    strip_re = re.compile(r'\t(?:\d+_%s_)?([^\t]*)%s\t' % (genstr, MASKTAG))

    # Keep this regex in sync with the file naming scheme used in split_chrs.
    subchr_re = re.compile(r'^(.*)_\+(\d+)$')
    with open(allpsl, 'wb') as allfh:
      for inp in psls:
        with open(inp, 'rb') as pfh:
          for line in pfh:

            # We need to rewrite the chrnames here. Also remove the
            # trailing newline so it doesn't confuse the processing below.
            newline = strip_re.sub(repl, line).rstrip('\n')

            # Parse out sub-chromosome coordinates from filenames and
            # fix coords appropriately. This is heavily dependent on
            # the PSL file following the specification.
            fields = newline.split("\t")
            if len(fields) > 1:

              # Sort out the query positions.
              chrA_match = subchr_re.match(fields[9])
              if chrA_match:
                fields[9] = chrA_match.group(1)
                basecoord = int(chrA_match.group(2))
                for fnum in (11,12):
                  fields[fnum] = str(int(fields[fnum]) + basecoord)
                fields[19] = ','.join([ str(int(x) + basecoord)
                                        for x in fields[19].split(',')
                                        if x != '']) + ','
                fields[10] = from_sizes[ fields[9] ]

              # Sort out the target positions.
              chrB_match = subchr_re.match(fields[13])
              if chrB_match:
                fields[13] = chrB_match.group(1)
                basecoord  = int(chrB_match.group(2))
                for fnum in (15,16):
                  fields[fnum] = str(int(fields[fnum]) + basecoord)
                fields[20] = ','.join([ str(int(x) + basecoord)
                                        for x in fields[20].split(',')
                                        if x != '']) + ','
                fields[14] = to_sizes[ fields[13] ]

              # Quick check on our output. This is essentially cribbed
              # from the pslToBed code.
              if (    int(fields[11]) >= int(fields[12])
                   or int(fields[12]) >  int(fields[10])
                   or int(fields[15]) >= int(fields[16])
                   or int(fields[16]) >  int(fields[14]) ):
                raise StandardError(
                  ("Mangled PSL format output. Offending input line was in file %s:"
                   + "\n\n%s\n\nMunged to:\n%s\n\n") % (inp, line, "\t".join(fields)))

            newline = "\t".join(fields) + "\n"
            
            allfh.write(newline)
        os.unlink(inp) # Attempt to save some temp space

    # Swap target and source annotation, such that splitting on the
    # target actually splits on the query.
    swppsl = os.path.join(self.local_tempdir, 'all-swap.psl')
    cmd = [ 'pslSwap', allpsl, swppsl ]
    call_subprocess(cmd, tmpdir=self.local_tempdir, path=os.environ['PATH'])
    os.unlink(allpsl)

    # Split psl files by target chromosome.
    psldir = os.path.join(self.local_tempdir, 'psl/')
    os.mkdir(psldir)

    # Consider -lump option for scaffolds FIXME
    cmd = [ 'pslSplitOnTarget', swppsl, psldir ]
    call_subprocess(cmd, tmpdir=self.local_tempdir, path=os.environ['PATH'])
    target_psls = [ os.path.join(self.local_tempdir, psldir, x)
                    for x in os.listdir(psldir) ]
    self._tempfiles.extend(target_psls + [psldir])
    os.unlink(swppsl)

    return target_psls

  def get_chr_sizes(self, fasta):
    '''
    Runs faSize on a fasta file to generate chr size data.
    '''

    # We keep a cached because we'll be using this more than once.
    if fasta in self._chr_sizes:
      return self._chr_sizes[fasta]
    LOGGER.info("Calculating chr sizes for %s", fasta)
    sizefn = os.path.join(self.local_tempdir, filebasename(fasta) + '.sizes')
    cmd = 'faSize %s -detailed > %s' % (bash_quote(fasta), bash_quote(sizefn))
    call_subprocess(cmd, tmpdir=self.local_tempdir,
                    path=os.environ['PATH'], shell=True)
    self._tempfiles.append(sizefn)
    self._chr_sizes[fasta] = sizefn
    return sizefn

  def get_chr_sizes_dict(self, fasta):
    '''
    As for get_chr_sizes, but also parses the file and returns a dict
    for convenience.
    '''
    sizefn = self.get_chr_sizes(fasta)
    sizes = dict()
    with open(sizefn, 'r') as sizefh:
      for row in sizefh:
        (chrom, size) = [ x.strip() for x in row.split() ]
        sizes[chrom] = size
    return sizes

  def chain(self, lavs):
    '''
    Chains the lastz output .lav files together.
    '''
    # We keep the filtered chain file.
    gen_from = filebasename(self.from_genome)
    gen_to   = filebasename(self.to_genome)
    prechain = os.path.join(self.local_workdir,
                            '%s_vs_%s.pre.chain' % (gen_from, gen_to))
    if os.path.exists(prechain):
      LOGGER.warning("Prechain file already exists."
                     + " Assuming we can start from this point: %s", prechain)
      return prechain

    # Convert lavs to appropriately-organised psl files.
    psls = self.process_lavs_to_psl(lavs)

    # FIXME at some point we need to add these psls to self._tempfiles

    # Run the initial chaining.
    LOGGER.info("Running the initial chaining.")
    chaindir = os.path.join(self.local_tempdir, 'chain/')
    os.mkdir(chaindir)
    chains = []
    for psl in psls:
      chfn = os.path.join(chaindir, filebasename(psl) + '.chain')
      cmd = [ 'axtChain', '-psl', '-linearGap=%s' % self.linear_gap, psl,
              '-faQ', self.from_genome,
              '-faT', self.to_genome, chfn ]
      call_subprocess(cmd, tmpdir=self.local_tempdir, path=os.environ['PATH'])
      chains.append(chfn)
      self._tempfiles.append(chfn)
    self._tempfiles.append(chaindir)

    # Filter the chained alignments before returning.
    allchain = os.path.join(self.local_tempdir, 'all.chain')
    cmd = ('chainMergeSort -tempDir=%s %s > %s'
           % (bash_quote(self.local_tempdir),
              " ".join([ bash_quote(x) for x in chains ]),
              bash_quote(allchain)))
    call_subprocess(cmd, tmpdir=self.local_tempdir,
                    path=os.environ['PATH'], shell=True)
    self._tempfiles.append(allchain)

    from_sizes = self.get_chr_sizes(self.from_genome)
    to_sizes   = self.get_chr_sizes(self.to_genome)

    # Actually create the prechain file.
    cmd = [ 'chainPreNet', allchain, from_sizes, to_sizes, prechain ]
    call_subprocess(cmd, tmpdir=self.local_tempdir, path=os.environ['PATH'])
    return prechain

  def get_2bit(self, fasta):
    '''
    Simply generate a temporary 2bit file from the specified fasta
    file. Note the differences between this and convert_to_2bit. FIXME
    refactor so there's only one of these functions.
    '''
    outfn = os.path.join(self.local_tempdir, filebasename(fasta) + '.2bit')
    if os.path.exists(outfn):
      return outfn
    LOGGER.info("Generating 2bit file for %s", fasta)
    cmd = [ 'faToTwoBit', fasta, outfn ]
    call_subprocess(cmd, tmpdir=self.local_tempdir, path=os.environ['PATH'])
    self._tempfiles.append(outfn)
    return outfn

  def net(self, prechain):
    '''
    Create nets from the chained alignements and convert them to axt
    format. Also generate a liftOver file.
    '''
    from_sizes = self.get_chr_sizes(self.from_genome)
    to_sizes   = self.get_chr_sizes(self.to_genome)

    net = os.path.join(self.local_workdir, prechain + '.net')
    cmd = (('chainNet %s -minSpace=1 %s %s stdout /dev/null'
            + ' | netSyntenic stdin %s')
           % (bash_quote(prechain), bash_quote(from_sizes),
              bash_quote(to_sizes), bash_quote(net)))
    # This may fail for spurious reasons (e.g. absence of
    # /proc/self/stat on non-linux machines).
    try:
      LOGGER.info("Running chainNet and netSyntenic on prechain file.")
      call_subprocess(cmd, tmpdir=self.local_tempdir,
                      path=os.environ['PATH'], shell=True)
    except CalledProcessError, err:
      LOGGER.warning("chainNet or netSyntenic raised exception: %s", err)
    if not os.path.exists(net):
      raise StandardError(
        "chainNet/netSyntenic failed to create output net file %s" % net)

    axt = os.path.join(self.local_workdir,
                       "%s.%s.net.axt" % (filebasename(self.from_genome),
                                          filebasename(self.to_genome)))
    from_2bit = self.get_2bit(self.from_genome)
    to_2bit   = self.get_2bit(self.to_genome)
    LOGGER.info('Converting to axt format.')
    cmd = ('netToAxt %s %s %s %s stdout | axtSort stdin %s'
           % (bash_quote(net), bash_quote(prechain),
              bash_quote(from_2bit), bash_quote(to_2bit),
              bash_quote(axt)))
    call_subprocess(cmd, tmpdir=self.local_tempdir,
                    path=os.environ['PATH'], shell=True)

    # These are cheap to generate and store, but potentially very useful later.
    LOGGER.info('Creating liftOver file.')
    liftover = os.path.join(self.local_workdir, prechain + '.liftOver')
    cmd = ('netChainSubset',  net, prechain, liftover)
    call_subprocess(cmd, tmpdir=self.local_tempdir,
                    path=os.environ['PATH'])
    
    return axt

  def preprocess(self, fasta):
    '''
    Convenience method to prepare a single genome fasta file for
    alignment with lastz. Returns a list of masked 2bit files created
    in the designated working directory.
    '''
    LOGGER.info("Commencing preprocessing for file: %s", fasta)
    workdir = os.path.join(self.local_workdir, filebasename(fasta))
    if not os.path.exists(workdir):
      os.mkdir(workdir)

    # If there are any 2bit files in the target directory it's assumed
    # that the set is complete and no repeat preprocessing is
    # necessary.
    twobit_re = re.compile(r'\.2bit$')
    twobits   = [ os.path.join(workdir, x) for x in os.listdir(workdir)
                                         if twobit_re.search(x) ]
    if len(twobits) == 0:
      split  = self.split_chrs(fasta)
      self._tempfiles.extend(split)
      masked = [ self.mask_tandem_repeats(x) for x in split ]
      self._tempfiles.extend(masked)
      twobits  = [ self.convert_to_2bit( x, workdir ) for x in masked ]
    else:
      LOGGER.warning("Working directory %s contains 2bit files."
                     + " Assuming these are complete!", workdir)

    return twobits

  def generate_lav(self):
    '''
    Method runs the first part of the lastz alignment pipeline, up to
    and including the generation of .lav files by lastz on the
    cluster.
    '''
    lav_re = re.compile(r'\.lav$')
    lavfiles   = [ os.path.join(self.local_workdir, x)
                   for x in os.listdir(self.local_workdir)
                   if lav_re.search(x) ]

    if len(lavfiles) == 0 or self.resume:

      # Process Genome file 1 to per-chromosome masked 2bit files.
      from_2bits = self.preprocess(self.from_genome)

      # Process Genome file 2.
      to_2bits   = self.preprocess(self.to_genome)

      # Run the alignment.
      (jobs, lavfiles, cluster2bits) = \
          self.align(from_2bits, to_2bits,
                     [ os.path.basename(x) for x in lavfiles ])

      # Set up a notifier LSF job.
      # Set up listening socket, wait for that nc command to run.
      self.wait_on_cluster(jobs,
                           cleanup_cmd="rm %s" % " ".join(cluster2bits))
      self.retry_lavfile_scp()
      
    else:
      LOGGER.warning("Working directory %s contains lav files."
                     + " Will not re-run the alignment.", self.local_workdir)
    return lavfiles

  def retry_lavfile_scp(self):
    '''
    If any completed lavfiles failed to transfer, retry (just the
    once, in real time). Typically (I think) the lavfile scp step fails
    due to too many lastz jobs completing at once, so this slows the
    demand on localhost ssh port to a virtual trickle for the last set
    of transfers.
    '''
    cmd = r'find . -maxdepth 1 -name %d_\*.lav.done' % os.getpid()
    with self.runner.run_command(cmd) as query:
      orphans = [  os.path.splitext(x.strip())[0] for x in query ]
    for lav in orphans:
      outfile = re.sub(r'^\d+_', '', lav)
      if not os.path.exists(os.path.join(self.local_workdir, outfile)):
        LOGGER.info("Retrying scp of missing but completed file %s", lav)
        cmd = self.return_file_to_localhost(lav, outfile)

        # Also clean up the lav and lav.done files on the cluster.
        cmd = 'rm %s %s.done' % (lav, lav)
        self.runner.run_command(cmd)

  def run_pipeline(self):
    '''
    Run the pipeline as described on the UCSC genomewiki pages. Split
    the input genomes by chromosome, mask tandem repeats, convert to
    2bit format, run lastz (ideally via a cluster), chain and net the
    results.
    '''
    LOGGER.info("Starting alignment pipeline.")
    lavfiles = self.generate_lav()

    # Run chaining and netting.
    LOGGER.info("Commencing chaining.")
    prechain = self.chain(lavfiles)

    self.net(prechain)
    LOGGER.info("Analysis complete.")

  def __del__(self):
    '''
    Attempt to clean up our temporary files and directories.
    '''
    sys.stderr.write("Cleaning up temporary files.\n")
    for fname in self._tempfiles:
      if os.path.exists(fname):
        try:
          if os.path.isdir(fname):
            os.removedirs(fname)
          else:
            os.unlink(fname)
        except OSError, err:
          sys.stderr.write(
            "WARNING: Unable to delete temporary file or directory %s:\n%s\n"
            % (fname, err))
    try:
      os.removedirs(self.local_tempdir)
    except OSError, err:
      sys.stderr.write(
        "WARNING: Unable to delete temporary directory %s:\n%s\n"
        % (self.local_tempdir, err))

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(
    description='Runs the LASTZ alignment, chaining and netting'
    + ' pipeline from UCSC on two genomes')

  PARSER.add_argument('-f', '--from', dest='from_genome', type=str,
                      required=True, help='A genome to align.')

  PARSER.add_argument('-t', '--to', dest='to_genome', type=str,
                      required=True, help='A genome to align.')

  PARSER.add_argument('-d', '--dir', dest='dir', type=str, required=True,
                      help='The working directory in which intermediate'
                      + ' results are stored.')

  PARSER.add_argument('--hsp-thresh', dest='hsp_thresh',
                      type=int, default=3000,
                      help='The hsp-thresh argument to use with lastz'
                      + ' (default=3000; this is likely too low).')

  PARSER.add_argument('--length-limit', dest='length_limit',
                      type=int, default=None,
                      help='The chromosome length limit (in bases) above'
                      + ' which sequences will be split into sub-segments'
                      + ' for alignment (default=Unlimited).')

  PARSER.add_argument('--linear-gap', dest='linear_gap',
                      type=str, default='loose',
                      help='The linearGap option to use with axtChain'
                      + ' (loose (DEFAULT) for chicken->human, medium for'
                      + ' mouse->human; can be a filename specifying'
                      + ' linearGap parameters).')

  PARSER.add_argument('-i', '--ssh-key', dest='ssh_key', type=str,
                      help='The path to a (password-free) ssh key used by the'
                      + ' cluster account to transfer the completed alignment'
                      + ' files back to our local host.')

  PARSER.add_argument('--resume', dest='resume', action='store_true',
                      help='Indicates that the script should resume an'
                      + ' earlier run by filling in missing'
                      + ' lav files, rather than assuming that the list of lav'
                      + ' files in the working directory is complete.')

  PARSER.add_argument('--throttle', dest='throttle', type=int, default=0,
                      help='The number of seconds to wait between submitting'
                      + ' cluster jobs. This can be used to give other cluster'
                      + ' users a chance to get jobs submitted, and to maybe'
                      + ' reduce the total number of concurrent jobs running.')

  PARSER.add_argument('--memsize', dest='memsize', type=int, default=5,
                      help='The amount of memory, in megabytes, to request on'
                      + ' the cluster (default=5).')

  PARSER.add_argument('--timelimit', dest='time_limit', type=int, default=48,
                      help='The length of time, in hours, to allow cluster'
                      + ' jobs to run for (default=48 hours).')

  PARSER.add_argument('--dryrun', dest='dryrun', action='store_true',
                      help='Run the script in dry-run mode. Currently this'
                      + ' simply figures out how many cluster jobs would be run'
                      + ' for a given set of genomes and length_limit parameter.')

  ARGS = PARSER.parse_args()

  ALIGNER = LastzAligner(from_genome   = ARGS.from_genome,
                         to_genome     = ARGS.to_genome,
                         local_workdir = ARGS.dir,
                         hsp_thresh    = ARGS.hsp_thresh,
                         ssh_key       = ARGS.ssh_key,
                         length_limit  = ARGS.length_limit,
                         linear_gap    = ARGS.linear_gap,
                         throttle      = ARGS.throttle,
                         memsize       = ARGS.memsize,
                         time_limit    = ARGS.time_limit,
                         resume        = ARGS.resume)

  if ARGS.dryrun:

    # Keep this simple for now; do not embed confusing (and
    # hard-to-maintain) debugging code throughout the class!
    FROM = len(ALIGNER.split_chrs(ARGS.from_genome, dryrun=True))
    TO   = len(ALIGNER.split_chrs(ARGS.to_genome,   dryrun=True))
    print("Requested settings would run %d cluster jobs: %d (from) x %d (to)."
          % (FROM*TO, FROM, TO))
  else:

    # Run the actual pipeline.
    ALIGNER.run_pipeline()
