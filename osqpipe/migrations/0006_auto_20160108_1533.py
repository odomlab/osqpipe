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

from django.contrib.auth.models import User
from rest_framework.authtoken.models import Token

'''
Creates initial token set for all users.
'''

def create_user_tokens(apps, schema_editor):

    for user in User.objects.all():
        Token.objects.get_or_create(user=user)

class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0005_auto_20151218_1551'),
        ('authtoken', '0001_initial'),
    ]

    operations = [
        migrations.RunPython(create_user_tokens),
    ]
