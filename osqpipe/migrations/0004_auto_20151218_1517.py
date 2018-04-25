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
        ('osqpipe', '0003_auto_20151218_1504'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='genome',
            options={'ordering': ['code']},
        ),
        migrations.RemoveField(
            model_name='genome',
            name='common_name',
        ),
        migrations.RemoveField(
            model_name='genome',
            name='scientific_name',
        ),
        migrations.RemoveField(
            model_name='genome',
            name='taxonomy_id',
        ),
        migrations.AlterField(
            model_name='genome',
            name='species',
            field=models.ForeignKey(to='osqpipe.Species', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AlterField(
            model_name='species',
            name='accession',
            field=models.CharField(unique=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='species',
            name='scientific_name',
            field=models.CharField(unique=True, max_length=255, db_column=b'sciname'),
        ),
    ]
