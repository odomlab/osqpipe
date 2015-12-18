# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0004_auto_20151218_1517'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='characteristic',
            unique_together=set([('category', 'value')]),
        ),
    ]
