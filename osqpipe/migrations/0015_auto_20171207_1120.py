# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0014_restrictome'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='restrictome',
            name='genome',
        ),
        migrations.DeleteModel(
            name='Restrictome',
        ),
    ]
