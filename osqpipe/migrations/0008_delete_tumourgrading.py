# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0007_auto_20160302_1045'),
    ]

    operations = [
        migrations.DeleteModel(
            name='TumourGrading',
        ),
    ]
