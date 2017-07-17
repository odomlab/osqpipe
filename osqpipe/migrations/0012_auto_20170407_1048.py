# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0011_auto_20170316_1327'),
    ]

    operations = [
        migrations.CreateModel(
            name='Condition',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=255)),
                ('description', models.TextField(null=True, blank=True)),
            ],
            options={
                'ordering': ['name'],
                'db_table': 'condition',
            },
        ),
        migrations.AddField(
            model_name='library',
            name='condition',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.Condition', null=True),
        ),
    ]
