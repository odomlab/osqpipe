# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0012_auto_20170407_1048'),
    ]

    operations = [
        migrations.AddField(
            model_name='project',
            name='is_frozen',
            field=models.BooleanField(default=False),
        ),
    ]
