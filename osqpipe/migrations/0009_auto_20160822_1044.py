# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0008_delete_tumourgrading'),
    ]

    operations = [
        migrations.AddField(
            model_name='library',
            name='platecode',
            field=models.CharField(max_length=32, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='library',
            name='platecol',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='library',
            name='platerow',
            field=models.CharField(blank=True, max_length=1, null=True, choices=[(b'A', b'A'), (b'B', b'B'), (b'C', b'C'), (b'D', b'D'), (b'E', b'E'), (b'F', b'F'), (b'G', b'G'), (b'H', b'H')]),
        ),
    ]
