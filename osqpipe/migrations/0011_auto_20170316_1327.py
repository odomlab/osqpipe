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


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0010_alignmentqc_alnqcfile'),
    ]

    operations = [
        migrations.AddField(
            model_name='library',
            name='external_records',
            field=models.ManyToManyField(related_name='libraries', db_table=b'library_external_record', to='osqpipe.ExternalRecord'),
        ),
        migrations.AddField(
            model_name='sample',
            name='external_records',
            field=models.ManyToManyField(related_name='samples', db_table=b'sample_external_record', to='osqpipe.ExternalRecord'),
        ),
    ]
