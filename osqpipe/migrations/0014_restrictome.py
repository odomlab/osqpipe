# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0013_project_is_frozen'),
    ]

    operations = [
        migrations.CreateModel(
            name='Restrictome',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('enzyme', models.CharField(max_length=128)),
                ('sequence', models.CharField(max_length=128, null=True, blank=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('date', models.DateField(auto_now_add=True)),
                ('version', models.CharField(max_length=255, null=True, blank=True)),
                ('genome', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='osqpipe.Genome', help_text=b'The genome against which the restrictome was generated.')),
            ],
            options={
                'ordering': ['enzyme', 'genome'],
                'db_table': 'restrictome',
            },
        ),
    ]
