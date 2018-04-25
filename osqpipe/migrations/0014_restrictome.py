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
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0013_project_is_frozen'),
    ]

    operations = [
        migrations.CreateModel(
            name='Restrictome',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('enzyme', models.CharField(max_length=128)),
                ('sequence', models.CharField(max_length=128, null=True, blank=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('date', models.DateField(auto_now_add=True)),
                ('version', models.CharField(max_length=255, null=True, blank=True)),
                ('genome', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='osqpipe.Genome', help_text=b'The genome against which the restrictome was generated.')),
            ],
            options={
                'ordering': ['enzyme', 'genome'],
                'db_table': 'restrictome',
            },
        ),
    ]
