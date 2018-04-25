=======
OSQPipe
=======

This package, "osqpipe", contains the code used to run the Odom lab
sequencing data analysis pipeline.

The package has been released as-is, and comes with no promise of any
support from its authors. The code was developed specifically to
address the needs of the Odom Lab during its residence at the CRUK
Cambridge Institute, and as such contains many optimisations which may
not be fit for purpose in other computational infrastructure
environments. Please see the LICENSE.txt file for further disclaimers.

Quick start
-----------

These instructions assume that you have set up a site-specific Django
project in which you are going to install the osqpipe app. Plese
refer to the Django documentation for help with this
(https://docs.djangoproject.com - see in particular the
tutorial). Once you have set up your project the following
instructions will hopefully make more sense.

1. Add "osqpipe" to your INSTALLED_APPS setting like this::

    INSTALLED_APPS = (
        ...
        'osqpipe',
    )

2. Include the osqpipe repository URLconf in your project urls.py like this::

    url(r'^repository/', include('osqpipe.urls')),

3. Run `python manage.py migrate` to create the osqpipe models.

4. Start the development server and visit http://127.0.0.1:8000/admin/
   to set up users etc. (you'll need the Admin app enabled).

5. Visit http://127.0.0.1:8000/repository/ to access the osqpipe repository database.

External Prerequisites
----------------------

Firstly, the following non-core python modules are required for this
package to function:

   * Django
   * django-dbarray
   * gnuplot-py
   * requests
   * fuzzy

Optional modules:

   * xlrd
   * BeautifulSoup4
   * lxml
   * pexpect
   * pysam (reallocateReads)
   * histogram, sequence, fasta, fastq, iset, some other CRI modules by Gord (?)

Secondly, this pipeline relies on many external programs to perform
its tasks. Most of the data management tasks are performed on the
primary host. Sequence alignment against reference genomes, and some
accompanying tasks, are run either on an LSF-based cluster or on a
single host with many cores available for multithreaded operation. The
programs needed on each server are as follows:

1. Primary host (config option: hostpath)::

   * python
   * samtools
   * fastqc
   * ssh
   * scp

   The following are part of the kent src utilities::

      * bedGraphToBigWig
      * fetchChromSizes

   The following are part of the Kraken suite::

      * reaper
      * tally

   The following are maintained by CRUK-CI::
   
      * demuxIllumina
      * export2fastq
      * makeWiggle
      * solexa2phred
      * summarizeFile

      * reallocateReads
      * trimFastq
      
      * bam2fq (thin wrapper for picard-tools; FIXME use those directly)

2. LSF cluster (config option: clusterpath)::

   * python
   * bsub
   * samtools
   * bwa
   * split
   * scp
   * gzip
   * sh (i.e., bash)

3. Multicore alignment host (config option: althostpath)::

   * python
   * samtools
   * bwa
   * ssh
   * scp
   * gzip
   * rm  (FIXME review this also??)
   * sh (i.e., bash), nice, nohup

Note that the presence of a bash shell is assumed on both the LSF
cluster and the multicore alignment host.

Configuration
-------------

Database configuration is handled within your site-specific project
settings.py file as described in the Django documentation. For
osqpipe-specific settings (hostnames, various directories) you will
need to edit the file "osqpipe_config.xml". The pipeline will look
for this file in the following places, in order: current working
directory, user's home directory, /etc, and the directory pointed to
by the $OSQPIPE_CONFDIR environmental variable.

FIXME needs an explanation of the various config settings, either
here, or in "hint" attributes in the XML config itself.

Credits
-------

The code in this package is the product of a joint sequencing effort
by Gord Brown, Tim Rayner and Margus Lukk during the years of Odom Lab
operation at the CRUK Cambridge Institute. The sequencing pipeline
code in this package was originally developed by Gord Brown working
with the Odom and Carroll labs. Additional features were implemented
by Margus Lukk and Tim Rayner. The codebase was refactored and
migrated to use the Django framework by Tim Rayner.
