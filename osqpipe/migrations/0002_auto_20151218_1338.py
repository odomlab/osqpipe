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

from django.db import models, migrations
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Characteristic',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('category', models.CharField(max_length=32)),
                ('value', models.CharField(max_length=32)),
            ],
            options={
                'ordering': ['category', 'value'],
                'db_table': 'characteristic',
            },
        ),
        migrations.CreateModel(
            name='SizeUnit',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=32)),
                ('description', models.CharField(max_length=128)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'size_unit',
            },
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('scientific_name', models.CharField(max_length=255, db_column=b'sciname')),
                ('common_name', models.CharField(max_length=255, null=True, db_column=b'commonname', blank=True)),
                ('accession', models.CharField(max_length=32)),
            ],
            options={
                'ordering': ['scientific_name'],
                'db_table': 'species',
                'verbose_name_plural': 'species',
            },
        ),
        migrations.RemoveField(
            model_name='sample',
            name='tumour_grading',
        ),
        migrations.AddField(
            model_name='sample',
            name='size',
            field=models.DecimalField(null=True, max_digits=5, decimal_places=2, blank=True),
        ),
        migrations.AddField(
            model_name='genome',
            name='species',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Species', null=True),
        ),
        migrations.AddField(
            model_name='sample',
            name='characteristics',
            field=models.ManyToManyField(related_name='samples', db_table=b'sample_characteristic', to='osqpipe.Characteristic'),
        ),
        migrations.AddField(
            model_name='sample',
            name='size_unit',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.SizeUnit', null=True),
        ),
        migrations.AddField(
            model_name='source',
            name='species',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Species', null=True),
        ),
    ]
