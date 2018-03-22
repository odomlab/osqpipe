# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0016_restrictome'),
    ]

    operations = [
        migrations.AddField(
            model_name='library',
            name='release_worthy',
            field=models.BooleanField(default=False),
        ),
    ]
