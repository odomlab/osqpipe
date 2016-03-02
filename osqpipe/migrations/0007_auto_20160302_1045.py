# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0006_auto_20160108_1533'),
    ]

    operations = [
        migrations.AlterField(
            model_name='sample',
            name='name',
            field=models.CharField(max_length=128),
        ),
        migrations.AlterField(
            model_name='source',
            name='name',
            field=models.CharField(unique=True, max_length=128),
        ),
    ]
