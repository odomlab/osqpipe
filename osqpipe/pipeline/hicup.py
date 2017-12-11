#!/usr/bin/env python
#
# $Id$
#
import os
import sys
from subprocess import Popen, PIPE
import shutil

from osqpipe.models import Restrictome, Genome, Program
import django
django.setup()

from osqutil.cluster import ClusterJobSubmitter
from osqutil.utilities import write_to_remote_file, transfer_file
from osqutil.config import Config


class HiCUP(object):
    
    def __init__(self, fq1, genome=None, enzyme=None, fq2=None):
        
        self.fq1 = fq1
        self.fq2 = fq2
        self.genome = genome
        self.genome_index = None
        self.enzyme = enzyme
        self.restriction_file = None

        self._check_file(fq1)
        self._check_file(fq2)

        self.alignment_program = 'bowtie2'
        
        self.conf = Config()

        self.hicup_conf_fname = os.path.join(self.conf.clusterworkdir, os.path.basename(self.fq1) + "_hicup.conf")
        self.hicup_output_dir = os.path.join(self.conf.clusterworkdir, os.path.basename(self.fq1) + "_hicup")

        # Get genome_file
        if self.genome is not None:
            self.genome_index = self._genome_index_path(genome)
            if enzyme is not None:
                self.restriction_file = self._restriction_file_path(genome, enzyme)
                
    def add_genome_file(self, genome_file):
        self._check_file(genome_file)
        self.genome_file = genome_file

    def add_restriction_file(self, restriction_file):
        self._check_file(restriction_file)
        self.restriction_file = restriction_file
        
    def _check_file(self, fname):

        if fname is not None:
            if not os.path.isfile(fname):
                LOGGER.error("%s not found or accessible! Exiting!", fname)
                sys.exit(1)

    def _genome_index_path(self, genome_code):
        '''Alternatively to the code below, one could use genome_fasta_path from .bwa_runner'''
        
        try:
            genome_obj = Genome.objects.get(code=genome)
        except Genome.DoesNotExist:
            LOGGER.error("Genome '%s' not found!", genome)
            sys.exit(1)
        try:
            prog_obj = Program.objects.get(program=self.alignment_program, current=True)
        except Program.DoesNotExist:
            LOGGER.error("Program '%s' not found!", self.alignment_program)
            sys.exit(1)

        sciname = genome_obj.species.scientific_name
        sciname = sciname.replace(" ", "_")
        sciname = sciname.lower()
        if self.alignment_program == "bowtie2":
            indexdir = "bowtie-%s" % prog_obj.version
        else:
            indexdir = "%s-%s" % (self.alignment_program, prog_obj.version)

        index_path = os.path.join(self.conf.clustergenomedir, sciname,
                             genome.code, indexdir, genome.code)

        return index_path

    def _restriction_file_path(genome, enzyme):

        try:
            genome_obj = Genome.objects.get(code=genome)
        except Genome.DoesNotExist:
            LOGGER.error("Genome '%s' not found!", genome)
            sys.exit(1)
        try:
            prog_obj = Program.objects.get(program='hicup', current=True)
        except Program.DoesNotExist:
            LOGGER.error("Program '%s' not found!", self.alignment_program)
            sys.exit(1)        
        try:
            r = Restrictome.objects.get(genome=genome_obj, enzyme=enzyme, program=prog_obj)
        except Restrictome.DoesNotExist:
            LOGGER.error("No restrictome for genome '%s' and enzyme '%s' for HiCUP", self.genome, self.enzyme)
            sys.exit(1)
        sciname = genome_obj.species.scientific_name
        sciname = sciname.replace(" ", "_")
        sciname = sciname.lower()
        indexdir = "%s-%s" % (self.alignment_program, prog_obj.version)    
        # E.g. /scratchb/user/fnc-odompipe/genomes/mus_musculus/mm10/hicup-0.5.10/Digest_mm10_HindIII_None_16-49-06_29-11-2017.txt
        rfile_path = os.path.join(self.conf.clustergenomedir, sciname, genome.code, indexdir, r.filename)

        return rfile_path

    def write_hicup_config(self):


        threads = int(int(self.conf.num_threads)/2)
        if self.fq2 is None:
            threads = self.conf.num_threads
        
        conf_txt =  "Outdir: %s\n" % self.hicup_output_dir
        conf_txt += "Threads: %s\n" % threads
        conf_txt += "Quiet:0\n"
        conf_txt += "Keep:0\n"
        conf_txt += "Zip:1\n"
        conf_txt += "Bowtie2: %s\n" % self.alignment_program
        conf_txt += "Index: %s\n" % self.genome_index
        conf_txt += "Digest: %s\n" % self.restriction_file
        conf_txt += "Format: Sanger\n"
        conf_txt += "Longest: 800\n"
        conf_txt += "Shortest: 150\n"
        conf_txt += "%s" % self.fq1
        if self.fq2 is not None:
            conf_txt += "%s\n" % self.fq2

        try:
            sshkey = self.conf.clustersshkey
        except AttributeError, _err:
            sshkey = None
            write_to_remote_file(conf_txt, self.hicup_conf_fname, self.conf.clusteruser,
                                 self.conf.cluster, append=False, sshkey=sshkey)
            
        return conf_fname

    def run_hicup(self, hicup_config):
        
        # Copy files
        # NB! There is vulnerability in below as we asssume input file follows odom lab convention
        code = self.fq1.split('_')[0]
        destination = "%s@%s:%s" % (self.conf.user, self.conf.cluster, self.conf.clusterworkdir)
        transfer_file(self.fq1, destination)
        if self.fq2 is not None:
            transfer_file(self.fq1, destination)

        # Send bsub for running hicup with correct number of thread request
        submitter = ClusterJobSubmitter()
        cmd = "cd %s && hicup --config %s" % (self.conf.clusterworkdir,self.hicup_conf_fname)
        jobid = submitter.submit_command(cmd=cmd, mem=self.conf.clustermem, auto_requeue=False, threads=self.conf.num_threads)
        
        #
        cmd = "cs_run_hicup_postprocess.py --fq1" % (self.conf.clusterworkdir,self.fq1)
        jobid = submitter.submit_command(cmd=cmd, mem=self.conf.clustermem, auto_requeue=False, threads=self.conf.num_threads, depend_jobs=[jobid])

        
    def postprocess_hicup(self):
        '''Post-processes hicup results. NB! The function is expected to run in cluster.'''
        
        # Find html report
        report_file = None
        if not os.path.isdir(self.hicup_output_dir):
            LOGGER.error("No report dir found! Expected %s.", self.hicup_output_dir)
            sys.exit(1)
        for f in os.listdir(self.hicup_output_dir):
            if f.endswith('html'):
                report_file = f
                break
        if report_file is None:
            LOGGER.error("No html report found in %s.", self.hicup_output_dir)
            sys.exit(1)

        # Copy report to repository
        # NB! There is vulnerability in below as we asssume input file follows odom lab convention
        code = self.fq1.split('_')[0]
        destination = "%s@%s:%s/%s/" % (self.conf.user, self.conf.datahost, self.conf.repositorydir, code)
        transfer_file(os.path.join(self.hicup_output_dir, report_file), destination)
    
        # Register report in repository
        md5 = checksum_file(report_file, unzip=False)
        cmd = "cs_addFile.py --qcfile --program_name hicup -M %s %s" % (os.path.join(destination, f), md5)

        subproc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        (stdout, stderr) = subproc.communicate()
        retcode = subproc.wait()
        if stdout is not None:
            LOGGER.info("STDOUT:")
            LOGGER.info(sys.stdout.write(stdout))
        if stderr is not None:
            LOGGER.error("STDERR:")
            LOGGER.error(sys.stderr.write(stderr))
        if retcode != 0:
            LOGGER.error("Failed to execute '%s'\n\n" % cmd)
            sys.exit(1)
        
        # Remove report dir
        shutil.rmtree(self.hicup_output_dir)
