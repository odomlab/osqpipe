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

#
# Created by lukk01 on 18/11/2016 for construction of BioSamples archive submission spreadsheets for MCE project.
#
import os
import sys
import time
import sqlite3

from osqpipe.models import Library, SourceTreatment, Characteristic
import django
django.setup()

from osqutil.setup_logs import configure_logging
from logging import INFO
LOGGER = configure_logging(level=INFO)

class TermSource(object):
    def __init__(self, fname):

        '''Reads term information from fname to sql table in memory. Serves info about different columns associated with the term'''

        # Filename for tab-delimited term source spreadsheet.
        # The spreadsheet is expected to contain header line followed by term information in following columns:
        # Term, Category, Term_Source_REF, Term_Source_ID, Term_Source_URI
        self.fname = fname
        self.con = None # handler for SQL lite table to hold content of name in memory
        self.cur = None # handler for interacting with SQL lite table with content from fname
        self.referenced_sources = dict() # dictionary of sources for which terms have been queried.

        if os.path.isfile(self.fname):        
            self._read_table_from_fname()
        else:
            LOGGER.error("File '%s' not found or unacessible. Exiting!" % self.fname)
            sys.exit(1)
        
    def _read_table_from_fname(self):

        '''Reads Term Source reference table from file to a SQL lite database in memory'''

        self.con = sqlite3.connect(':memory:')
        self.con.execute('''CREATE TABLE t (
        Term TEXT,
        Category TEXT,
        Term_Source_REF TEXT,
        Term_Source_ID TEXT,
        Term_Source_URI)''')
        fh = open(self.fname,'rb')
        header = True
        for line in fh:
            if header:
                header = False
                continue            
            line = line.rstrip('\n')
            cols = line.split('\t')
            self.con.execute("INSERT INTO t (Term, Category, Term_Source_REF, Term_Source_ID, Term_Source_URI) VALUES (?,?,?,?,?)",cols)
            self.con.commit()
        self.cur = self.con.cursor()
        
    def getSourceREF(self, term, category):

        '''Returns Term Source REF for a term + category pair'''

        if term == None or term =='':
            return ""
        self.cur.execute('SELECT Term_Source_REF, Term_Source_URI FROM t WHERE Term=\'%s\' AND Category=\'%s\'' % (term, category))
        res = self.cur.fetchone()
        source = ""
        if res == None:            
            LOGGER.error("ERROR! Record for term '%s' missing in %s. Exiting!" % (term, self.fname))
            sys.exit(1)
        else:
            source = res[0]
            # record that the source of question has been queried
            if res[0] != '':                
                self.referenced_sources[res[0]] = res[1]
#       print "TERM=\"%s\", CATEGORY=\"%s\", SOURCE=\"%s\"" % (term, category, source)        

        return source

    def getSourceID(self, term, category):

        '''Returns Term Source ID for term + category pair'''

        if term == None or term =='':
            return ""
        self.cur.execute('SELECT Term_Source_ID FROM t WHERE Term=\'%s\' AND Category=\'%s\'' % (term, category))
        id = self.cur.fetchone()
        if id == None:
            LOGGER.error("ERROR! Record for term '%s' missing in %s. Exiting!" % (term, self.fname))
            sys.exit(1)
        else:
            id = id[0]
#       print "TERM=\"%s\", CATEGORY=\"%s\", ID=\"%s\"" % (term, category, id)        

        return id

    def getReferencedSources(self):

        '''Returns dictionary of Term Source REFs and corresponding URIs that have been queried from the term source reference table.'''

        return self.referenced_sources

    # Following function is currenly not used!
    def printTermValues(self):

        '''Prints list of Terms in the Term Reference Table.'''

        for row in self.cur.execute('SELECT Term, count(Term) FROM t GROUP BY Term'):
            print row[0]

class BioSamples(object):

    '''A class for information extraction and formatting of all libraries part of a probject in BioSamples archive submission format.'''
    
    def __init__(self, project, termsource_file):

        self.project = project
        self.mci = ""
        self.scd = ""

        # Load reference Term Source information from table.
        self.ts = TermSource(termsource_file)    

        # Create SCD section of BioSamples spreadsheet
        self._create_SCD()
        # Create MSI section of BioSamples spreadsheet. Note that SCD section needs to be created first for Term Source information to be correctly populated.
        self._create_MSI()
        
    def getAnnotations(self,code):

        '''Returns annotations in dict for a given Odom lab library code'''
        
        annotations = dict()
    
        library = Library.objects.get(code=code)    
    
        annotations['Sample Name'] = library.sample.name
        annotations['Sample Accession'] = ''
        annotations['Sample Description'] = ''
        annotations['Organism'] = library.genome.species.scientific_name
        annotations['Sex'] = library.sample.source.sex.name
        annotations['strain'] = library.sample.source.strain.name
        annotations['project'] = self.project
        annotations['tissue_type'] = library.sample.tissue.name
        if library.sample.source.date_of_birth:
            annotations['date of birth'] = library.sample.source.date_of_birth.strftime("%Y-%m-%d %H:%M:%S.%f")
        else:
            LOGGER.error("ERROR date_of_birth for sample '%s' (%s) not found!" % (library.sample.name, code))
            annotations['date of birth'] = ""
        if library.sample.source.date_of_death:
            annotations['date of death'] = library.sample.source.date_of_death.strftime("%Y-%m-%d %H:%M:%S.%f")
        else:
            LOGGER.error("ERROR date_of_death for sample '%s' (%s) not found!" % (library.sample.name, code))
            annotations['date of death'] = ""

        # Check if sample has been treated
        try:
            so = SourceTreatment.objects.get(source=library.sample.source)
            annotations['treatment_agent'] = so.agent.name
            annotations['treatment_dose'] = so.dose
            if so.dose_unit.name == "ug/g body weight":
                annotations['Unit'] = "micrograms per gram of body weight"
            else:
                annotations['Unit'] = so.dose_unit.name
            annotations['treatment_date'] = so.date.strftime("%Y-%m-%d_%H:%M:%S.%f")
        except SourceTreatment.DoesNotExist:
            annotations['treatment_agent'] = ''
            annotations['treatment_dose'] = ''
            annotations['Unit'] = ''
            annotations['treatment_date'] = ''

        # Check if sample has diagnosis    
        try:
            annotations['diagnosis'] = Characteristic.objects.get(category__iexact='Diagnosis',samples=library.sample).value
        except Characteristic.DoesNotExist:
            annotations['diagnosis'] = ''

        # Check if sample has TumourGrade
        try:
            annotations['tumour_grading'] = Characteristic.objects.get(category__iexact='TumourGrade',samples=library.sample).value
        except Characteristic.DoesNotExist:
            annotations['tumour_grading'] = ''
        
        return annotations

    def _create_MSI(self):

        '''Creates MCI part of the BioSamples submission spreadsheet'''
        
        # Create MCI section with partial information ready for further editing in Excel.
        self.msi = '''[MSI]
Submission\tTitle\t\trequired Short title, 50 characters approx. Submission\tIdentifier\t\tMust be blank Assigned by BioSamples Database
Submission\tDescription\tx1\trequired Short description, one paragraph.
Submission\tVersion\t\toptional	Version of SampleTab specification (currently 1.2)
Submission\tReference Layer\t\tMust be blank If this submission is part of the reference layer, this will be "true". Otherwise it will be "false".
Submission\tUpdate Date\tx1\toptional Date this submission was last modified. Must be in a YYYY-MM-DD format.
Submission\tRelease Date\tx1\toptional Date to be made public on. If blank, it will be public immediately. Must be in a YYYY-MM-DD format.
Person Last Name\tFAMILYNAME\tRequired
Person Initials\t\tEither middle initial or first initial depending if a first name is present
Person First Name\tFIRSTNAME\t 	
Person Email\tEMAIL\t	
Person Role\tsubmitter\tShould be a child of role in EFO
Organization Name\tCancer Research UK Cambridge Institute\tRequired
Organization Address\tUniversity of Cambridge, Li Ka Shing Centre, Robinson Way, Cambridge, CB2 0RE, United Kingdom\tOne line, comma separated
Organization URI\thttp://www.cruk.cam.ac.uk/\tWeb site.
Organization Role\tinstitution\thould be a child of role in EFO
Publication PubMed ID\t\tValid PubMed ID, numeric only
Publication DOI\t\tValid Digital Object Identifier''' + "\n"

        # Add Term Source information for terms referenced in SCD section.
        tsn_cols = ["Term Source Name",]
        tsuri_cols = ["Term Source URI",]
        tsversion_cols = ["Term Source Version",]

        for key in self.ts.referenced_sources:
            tsn_cols.append(key)
            tsuri_cols.append(self.ts.referenced_sources[key])
            
        self.msi += "\t".join(tsn_cols) + "\n"
        self.msi += "\t".join(tsuri_cols) + "\n"
        self.msi += "\t".join(tsversion_cols) + "\n"

    def _create_SCD(self):
        '''Creates SCD part of the BioSamples submission spreadsheet'''
        
        sample_columns = ('Sample Name','Sample Accession','Sample Description')
        annotation_columns = ('Organism','Sex')
        characteristic_columns = ('strain','tissue_type','diagnosis','tumour_grading','treatment_agent')
        # "Characteristic[term]" which do not need to be followed by "Term Source REF" and "Term Source ID" columns
        characteristic_columns_nosource = ('treatment_dose','Unit','date of birth','date of death','treatment_date','project')


        # Print SCD header
        self.scd += "[SCD]\n"
        cols = []
        for col in sample_columns:
            cols.append(col)
        for col in annotation_columns:
            cols.append(col)
            cols.append("Term Source REF")
            cols.append("Term Source ID")
        for col in characteristic_columns:
            cols.append("Characteristic[" + col + "]")
            cols.append("Term Source REF")
            cols.append("Term Source ID")
        for col in characteristic_columns_nosource:
            cols.append("Characteristic[" + col + "]")
        self.scd += "\t".join(cols) + "\n"

        # Find libraries associated to projects
        project_libs = Library.objects.filter(projects__name=self.project)

        # For each library, print annotations
        for lib in project_libs:
            annotations = self.getAnnotations(lib.code)
            cols = []
            for col in sample_columns:
                cols.append(annotations[col])
            for col in annotation_columns:
                cols.append(annotations[col])
                cols.append(self.ts.getSourceREF(annotations[col],col))
                cols.append(self.ts.getSourceID(annotations[col],col))
            for col in characteristic_columns:
                cols.append(annotations[col])
                cols.append(self.ts.getSourceREF(annotations[col],col))
                cols.append(self.ts.getSourceID(annotations[col],col))
            for col in characteristic_columns_nosource:
                cols.append(annotations[col])
            self.scd += "\t".join(cols) + "\n"

    def write_table(self, fname='-'):
        '''Writes BioSamples submission spreadsheet to a file.'''

        if fname is '-':
            fh = sys.stdout
        else:
            fh = open(fname, 'wb')
        fh.write(self.msi)
        fh.write(self.scd)
            
if __name__ == '__main__':

  import argparse

  PARSER = argparse.ArgumentParser(
    description='Creates sample spreadsheet for a project in repository in a format for BioSamples archive submission.')

  PARSER.add_argument('-p', '--project', dest='project', type=str, required=True,
                      help='Name of the project in Odom lab repository. E.g. \'HCC_WGS\' or \'HCC_exome\'')

  PARSER.add_argument('-t', '--termsource_file', dest='termsource_file', type=str, default='MCE_project_TermSource_reference_table.txt',
                      help='A tab delimited file with term reference info in following columns: Term, Category, Term_Source_REF, Term_Source_ID, Term_Source_URI.')

  PARSER.add_argument('-o', '--outfile', dest='outfile', type=str, required=False, help='Output filename. Default: project.txt. If outfile is \'-\', file is written to STDOUT.')
  
  ARGS = PARSER.parse_args()

  # initiate BioSamples for a project
  bsamples = BioSamples(project=ARGS.project, termsource_file=ARGS.termsource_file)
  outfile = ARGS.project + '.txt'
  if ARGS.outfile:
      outfile = ARGS.outfile
  bsamples.write_table(fname=outfile)
