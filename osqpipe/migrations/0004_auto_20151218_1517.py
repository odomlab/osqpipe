# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0003_auto_20151218_1504'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='genome',
            options={'ordering': ['code']},
        ),
        migrations.RemoveField(
            model_name='genome',
            name='common_name',
        ),
        migrations.RemoveField(
            model_name='genome',
            name='scientific_name',
        ),
        migrations.RemoveField(
            model_name='genome',
            name='taxonomy_id',
        ),
        migrations.AlterField(
            model_name='genome',
            name='species',
            field=models.ForeignKey(to='osqpipe.Species', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AlterField(
            model_name='species',
            name='accession',
            field=models.CharField(unique=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='species',
            name='scientific_name',
            field=models.CharField(unique=True, max_length=255, db_column=b'sciname'),
        ),
    ]
