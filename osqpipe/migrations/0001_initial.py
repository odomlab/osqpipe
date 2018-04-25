# -*- coding: utf-8 -*-
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

from __future__ import unicode_literals

from django.db import migrations, models
import dbarray
import osqpipe.models
import django.db.models.deletion
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='LibraryExtra',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code_text_prefix', models.CharField(max_length=128, editable=False)),
                ('code_numeric_suffix', models.IntegerField(null=True, editable=False)),
            ],
            options={
                'db_table': 'library_extra',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Adapter',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=32)),
                ('sequence', models.CharField(max_length=32, null=True, blank=True)),
                ('protocol', models.CharField(max_length=32)),
            ],
            options={
                'ordering': ['code'],
                'db_table': 'adapter',
            },
        ),
        migrations.CreateModel(
            name='Alnfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('checksum', models.CharField(max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('archive_date', models.DateField(null=True, blank=True)),
            ],
            options={
                'ordering': ['filename'],
                'db_table': 'alnfile',
            },
        ),
        migrations.CreateModel(
            name='Antibody',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=255)),
                ('lot_number', models.CharField(default=b'unknown', max_length=64)),
                ('description', models.TextField(null=True, blank=True)),
            ],
            options={
                'ordering': ['name', 'lot_number'],
                'db_table': 'antibody',
                'verbose_name_plural': 'antibodies',
            },
        ),
        migrations.CreateModel(
            name='ArchiveLocation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=32)),
                ('root_path', models.CharField(unique=True, max_length=1024)),
                ('host', models.CharField(max_length=1024)),
                ('host_port', models.CharField(max_length=8)),
                ('host_path', models.CharField(max_length=1024)),
                ('host_user', models.CharField(max_length=128)),
                ('host_delete_timelag', models.IntegerField(null=True)),
            ],
            options={
                'db_table': 'archive_location',
            },
        ),
        migrations.CreateModel(
            name='DataProcess',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'db_table': 'data_process',
            },
        ),
        migrations.CreateModel(
            name='DataProvenance',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('rank_index', models.IntegerField()),
                ('parameters', models.CharField(max_length=255)),
            ],
            options={
                'db_table': 'data_provenance',
            },
        ),
        migrations.CreateModel(
            name='DoseUnit',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=32)),
                ('description', models.CharField(max_length=128)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'dose_unit',
            },
        ),
        migrations.CreateModel(
            name='ExternalRecord',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('accession', models.CharField(unique=True, max_length=32)),
                ('is_public', models.BooleanField(default=False)),
                ('release_date', models.DateField()),
            ],
            options={
                'ordering': ['accession'],
                'db_table': 'external_record',
            },
        ),
        migrations.CreateModel(
            name='ExternalRepository',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=32)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'external_repository',
            },
        ),
        migrations.CreateModel(
            name='Facility',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=10)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('description', models.TextField(null=True, blank=True)),
            ],
            options={
                'ordering': ['code'],
                'db_table': 'facility',
                'verbose_name_plural': 'facilities',
            },
        ),
        migrations.CreateModel(
            name='Factor',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('description', models.TextField(null=True, blank=True)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'factor',
            },
        ),
        migrations.CreateModel(
            name='Filetype',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=10)),
                ('name', models.CharField(unique=True, max_length=32)),
                ('description', models.TextField(null=True, blank=True)),
                ('suffix', models.CharField(unique=True, max_length=32, blank=True)),
                ('gzip', models.BooleanField(default=True)),
            ],
            options={
                'ordering': ['code'],
                'db_table': 'filetype',
            },
        ),
        migrations.CreateModel(
            name='Genome',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=32)),
                ('common_name', models.CharField(max_length=255, db_column=b'commonname')),
                ('scientific_name', models.CharField(max_length=255, db_column=b'sciname')),
                ('blastdb', models.CharField(max_length=255, null=True, blank=True)),
                ('fasta', models.CharField(max_length=255, null=True, blank=True)),
                ('fasta_md5sum', models.CharField(max_length=32, null=True, blank=True)),
                ('notes', models.TextField(null=True, blank=True)),
                ('version', models.CharField(max_length=255, null=True, blank=True)),
                ('url', models.CharField(max_length=256, null=True, blank=True)),
                ('taxonomy_id', models.IntegerField()),
            ],
            options={
                'ordering': ['common_name'],
                'db_table': 'genome',
            },
        ),
        migrations.CreateModel(
            name='HistologyImagefile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('checksum', models.CharField(max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('archive_date', models.DateField(null=True, blank=True)),
                ('block', models.CharField(max_length=32)),
                ('batch', models.CharField(max_length=32)),
                ('archive', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.ArchiveLocation', null=True)),
                ('filetype', models.ForeignKey(to='osqpipe.Filetype', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'ordering': ['filename'],
                'db_table': 'histology_imagefile',
            },
        ),
        migrations.CreateModel(
            name='Lane',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('flowcell', models.CharField(max_length=32)),
                ('rundate', models.DateField()),
                ('reads', models.IntegerField(null=True, blank=True)),
                ('passedpf', models.IntegerField(null=True, blank=True)),
                ('lanenum', models.IntegerField()),
                ('flowlane', models.IntegerField()),
                ('paired', models.BooleanField(default=False)),
                ('readlength', models.IntegerField(null=True, blank=True)),
                ('mapped', models.IntegerField(null=True, blank=True)),
                ('seqsamplepf', models.TextField(blank=True)),
                ('seqsamplebad', models.TextField(blank=True)),
                ('qualmeanpf', dbarray.FloatArrayField(null=True)),
                ('qualstdevpf', dbarray.FloatArrayField(null=True)),
                ('qualmean', dbarray.FloatArrayField(null=True)),
                ('qualstdev', dbarray.FloatArrayField(null=True)),
                ('summaryurl', models.CharField(max_length=1024, null=True, blank=True)),
                ('genomicssampleid', models.CharField(max_length=32, null=True, blank=True)),
                ('usersampleid', models.CharField(max_length=1024, null=True, blank=True)),
                ('notes', models.TextField(null=True, blank=True)),
                ('failed', models.BooleanField(default=False)),
                ('runnumber', models.CharField(max_length=255, null=True, blank=True)),
                ('external_records', models.ManyToManyField(related_name='lanes', db_table=b'lane_external_record', to='osqpipe.ExternalRecord')),
                ('facility', models.ForeignKey(to='osqpipe.Facility', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'ordering': ['library'],
                'db_table': 'lane',
            },
        ),
        migrations.CreateModel(
            name='Lanefile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('checksum', models.CharField(max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('archive_date', models.DateField(null=True, blank=True)),
                ('pipeline', models.CharField(default=b'chipseq', max_length=128)),
                ('archive', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.ArchiveLocation', null=True)),
                ('filetype', models.ForeignKey(to='osqpipe.Filetype', on_delete=django.db.models.deletion.PROTECT)),
                ('lane', models.ForeignKey(to='osqpipe.Lane', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'ordering': ['filename'],
                'db_table': 'lanefile',
            },
        ),
        migrations.CreateModel(
            name='Library',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=128, validators=[osqpipe.models.validate_library_code])),
                ('bad', models.BooleanField(default=False)),
                ('barcode', models.CharField(max_length=32, null=True, blank=True)),
                ('chipsample', models.CharField(max_length=255, null=True, blank=True)),
                ('paired', models.BooleanField(default=False)),
                ('comment', models.TextField(null=True, blank=True)),
                ('adapter', models.ForeignKey(related_name='libraries', on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Adapter', null=True)),
                ('adapter2', models.ForeignKey(related_name='libraries2', on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Adapter', null=True)),
                ('antibody', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Antibody', null=True)),
                ('factor', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Factor', null=True)),
                ('genome', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='osqpipe.Genome', help_text=b'The genome against which the sequence data from this library should be aligned.')),
            ],
            options={
                'ordering': ['code'],
                'db_table': 'library',
                'verbose_name_plural': 'libraries',
            },
        ),
        migrations.CreateModel(
            name='LibraryNameMap',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('libname', models.CharField(max_length=128)),
                ('limsname', models.CharField(unique=True, max_length=128)),
            ],
            options={
                'db_table': 'library_name_map',
            },
        ),
        migrations.CreateModel(
            name='Libtype',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=32)),
                ('name', models.CharField(unique=True, max_length=32)),
                ('description', models.CharField(max_length=255, null=True, blank=True)),
            ],
            options={
                'ordering': ['code'],
                'db_table': 'libtype',
            },
        ),
        migrations.CreateModel(
            name='Linkerset',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('fivep', models.CharField(max_length=255)),
                ('threep', models.CharField(max_length=255)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'linkerset',
            },
        ),
        migrations.CreateModel(
            name='Machine',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=32)),
                ('platform', models.CharField(max_length=32)),
                ('name', models.CharField(max_length=32)),
            ],
            options={
                'ordering': ['code'],
                'db_table': 'machine',
                'verbose_name_plural': 'machines',
            },
        ),
        migrations.CreateModel(
            name='MergedAlnfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('checksum', models.CharField(max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('archive_date', models.DateField(null=True, blank=True)),
                ('archive', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.ArchiveLocation', null=True)),
                ('filetype', models.ForeignKey(to='osqpipe.Filetype', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'ordering': ['filename'],
                'db_table': 'merged_alnfile',
            },
        ),
        migrations.CreateModel(
            name='Peakfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('checksum', models.CharField(max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('archive_date', models.DateField(null=True, blank=True)),
                ('archive', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.ArchiveLocation', null=True)),
                ('filetype', models.ForeignKey(to='osqpipe.Filetype', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'ordering': ['filename'],
                'db_table': 'peakfile',
            },
        ),
        migrations.CreateModel(
            name='Program',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('program', models.CharField(max_length=128)),
                ('version', models.CharField(max_length=128)),
                ('options', models.CharField(max_length=256, null=True, blank=True)),
                ('files', models.CharField(max_length=256, null=True, blank=True)),
                ('type', models.CharField(max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('current', models.BooleanField(default=True)),
            ],
            options={
                'ordering': ['program', 'version'],
                'db_table': 'program',
            },
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=32)),
                ('name', models.CharField(unique=True, max_length=32)),
                ('description', models.TextField(null=True, blank=True)),
                ('shortnames', models.BooleanField(default=True)),
                ('filtered', models.BooleanField(default=True)),
                ('lab', models.CharField(max_length=32)),
                ('people', models.ManyToManyField(to=settings.AUTH_USER_MODEL, db_table=b'project_users')),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'project',
            },
        ),
        migrations.CreateModel(
            name='QCfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('checksum', models.CharField(max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('archive_date', models.DateField(null=True, blank=True)),
                ('archive', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.ArchiveLocation', null=True)),
                ('filetype', models.ForeignKey(to='osqpipe.Filetype', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'ordering': ['filename'],
                'db_table': 'qcfile',
                'verbose_name': 'QC file',
            },
        ),
        migrations.CreateModel(
            name='QCValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=32)),
                ('value', models.CharField(max_length=32)),
            ],
            options={
                'db_table': 'qc_value',
                'verbose_name': 'QC Value',
            },
        ),
        migrations.CreateModel(
            name='Sample',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=64)),
                ('comment', models.TextField(null=True, blank=True)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'sample',
            },
        ),
        migrations.CreateModel(
            name='Sex',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=32)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'sex',
            },
        ),
        migrations.CreateModel(
            name='Source',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=64)),
                ('date_of_birth', models.DateField(null=True, blank=True)),
                ('date_of_death', models.DateField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
                ('father', models.ForeignKey(related_name='child_as_father', on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Source', null=True)),
                ('mother', models.ForeignKey(related_name='child_as_mother', on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Source', null=True)),
                ('sex', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Sex', null=True)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'source',
            },
        ),
        migrations.CreateModel(
            name='SourceTreatment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('date', models.DateField()),
                ('dose', models.CharField(max_length=128, null=True, blank=True)),
            ],
            options={
                'ordering': ['date', 'agent'],
                'db_table': 'source_treatment',
            },
        ),
        migrations.CreateModel(
            name='Status',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=32)),
                ('description', models.CharField(max_length=128)),
                ('colour', models.CharField(default=b'#FFFFFF', max_length=32)),
                ('lanerelevant', models.BooleanField(default=False)),
                ('libraryrelevant', models.BooleanField(default=False)),
                ('sortcode', models.IntegerField(default=0)),
                ('authority', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Facility', null=True)),
            ],
            options={
                'ordering': ['code'],
                'db_table': 'status',
                'verbose_name_plural': 'status flags',
            },
        ),
        migrations.CreateModel(
            name='Strain',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('description', models.TextField(null=True, blank=True)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'strain',
            },
        ),
        migrations.CreateModel(
            name='Tissue',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('description', models.TextField(null=True, blank=True)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'tissue',
            },
        ),
        migrations.CreateModel(
            name='TreatmentAgent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=64)),
                ('description', models.CharField(max_length=256)),
                ('accession', models.CharField(max_length=32)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'treatment_agent',
            },
        ),
        migrations.CreateModel(
            name='TumourGrading',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=64)),
                ('description', models.CharField(max_length=256)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'tumour_grading',
            },
        ),
        migrations.CreateModel(
            name='Alignment',
            fields=[
                ('dataprocess_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='osqpipe.DataProcess')),
                ('total_reads', models.IntegerField()),
                ('mapped', models.IntegerField()),
                ('munique', models.IntegerField(null=True, blank=True)),
                ('headtrim', models.IntegerField(default=0, blank=True)),
                ('tailtrim', models.IntegerField(default=0, blank=True)),
                ('genome', models.ForeignKey(to='osqpipe.Genome', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'db_table': 'alignment',
            },
            bases=('osqpipe.dataprocess',),
        ),
        migrations.CreateModel(
            name='LaneQC',
            fields=[
                ('dataprocess_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='osqpipe.DataProcess')),
            ],
            options={
                'db_table': 'lane_qc',
                'verbose_name': 'Lane QC',
            },
            bases=('osqpipe.dataprocess',),
        ),
        migrations.CreateModel(
            name='MergedAlignment',
            fields=[
                ('dataprocess_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='osqpipe.DataProcess')),
                ('alignments', models.ManyToManyField(to='osqpipe.Alignment')),
            ],
            options={
                'db_table': 'merged_alignment',
            },
            bases=('osqpipe.dataprocess',),
        ),
        migrations.CreateModel(
            name='Peakcalling',
            fields=[
                ('dataprocess_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='osqpipe.DataProcess')),
                ('code', models.CharField(unique=True, max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('factor_align', models.ForeignKey(related_name='+', db_column=b'flane_id', on_delete=django.db.models.deletion.PROTECT, to='osqpipe.Alignment')),
                ('input_align', models.ForeignKey(related_name='+', db_column=b'ilane_id', on_delete=django.db.models.deletion.PROTECT, to='osqpipe.Alignment')),
            ],
            options={
                'db_table': 'peakcalls',
            },
            bases=('osqpipe.dataprocess',),
        ),
        migrations.AddField(
            model_name='sourcetreatment',
            name='agent',
            field=models.ForeignKey(to='osqpipe.TreatmentAgent', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='sourcetreatment',
            name='dose_unit',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.DoseUnit', null=True),
        ),
        migrations.AddField(
            model_name='sourcetreatment',
            name='source',
            field=models.ForeignKey(to='osqpipe.Source', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='source',
            name='strain',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Strain', null=True),
        ),
        migrations.AddField(
            model_name='sample',
            name='source',
            field=models.ForeignKey(to='osqpipe.Source', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='sample',
            name='tissue',
            field=models.ForeignKey(to='osqpipe.Tissue', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='sample',
            name='tumour_grading',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.TumourGrading', null=True),
        ),
        migrations.AlterUniqueTogether(
            name='program',
            unique_together=set([('program', 'version')]),
        ),
        migrations.AddField(
            model_name='library',
            name='libtype',
            field=models.ForeignKey(to='osqpipe.Libtype', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='library',
            name='linkerset',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Linkerset', null=True),
        ),
        migrations.AddField(
            model_name='library',
            name='projects',
            field=models.ManyToManyField(related_name='libraries', db_table=b'library_project', to='osqpipe.Project'),
        ),
        migrations.AddField(
            model_name='library',
            name='sample',
            field=models.ForeignKey(to='osqpipe.Sample', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='lane',
            name='library',
            field=models.ForeignKey(to='osqpipe.Library', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='lane',
            name='machine',
            field=models.ForeignKey(to='osqpipe.Machine', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='lane',
            name='status',
            field=models.ForeignKey(to='osqpipe.Status', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='histologyimagefile',
            name='sample',
            field=models.ForeignKey(to='osqpipe.Sample', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='externalrecord',
            name='repository',
            field=models.ForeignKey(to='osqpipe.ExternalRepository', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='dataprovenance',
            name='data_process',
            field=models.ForeignKey(related_name='provenance', to='osqpipe.DataProcess'),
        ),
        migrations.AddField(
            model_name='dataprovenance',
            name='program',
            field=models.ForeignKey(to='osqpipe.Program', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AlterUniqueTogether(
            name='antibody',
            unique_together=set([('name', 'lot_number')]),
        ),
        migrations.AddField(
            model_name='alnfile',
            name='archive',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.ArchiveLocation', null=True),
        ),
        migrations.AddField(
            model_name='alnfile',
            name='filetype',
            field=models.ForeignKey(to='osqpipe.Filetype', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AlterUniqueTogether(
            name='sourcetreatment',
            unique_together=set([('source', 'date', 'agent')]),
        ),
        migrations.AlterUniqueTogether(
            name='sample',
            unique_together=set([('name', 'tissue')]),
        ),
        migrations.AddField(
            model_name='qcvalue',
            name='laneqc',
            field=models.ForeignKey(related_name='qc_values', to='osqpipe.LaneQC'),
        ),
        migrations.AddField(
            model_name='qcfile',
            name='laneqc',
            field=models.ForeignKey(to='osqpipe.LaneQC', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='peakfile',
            name='peakcalling',
            field=models.ForeignKey(db_column=b'peakcalls_id', on_delete=django.db.models.deletion.PROTECT, to='osqpipe.Peakcalling'),
        ),
        migrations.AddField(
            model_name='mergedalnfile',
            name='alignment',
            field=models.ForeignKey(to='osqpipe.MergedAlignment', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='laneqc',
            name='lane',
            field=models.ForeignKey(to='osqpipe.Lane', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AlterUniqueTogether(
            name='lane',
            unique_together=set([('library', 'lanenum', 'facility')]),
        ),
        migrations.AlterUniqueTogether(
            name='dataprovenance',
            unique_together=set([('data_process', 'rank_index')]),
        ),
        migrations.AddField(
            model_name='alnfile',
            name='alignment',
            field=models.ForeignKey(to='osqpipe.Alignment', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='alignment',
            name='lane',
            field=models.ForeignKey(to='osqpipe.Lane', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AlterUniqueTogether(
            name='qcvalue',
            unique_together=set([('laneqc', 'name')]),
        ),
    ]
