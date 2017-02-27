# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0009_auto_20160822_1044'),
    ]

    operations = [
        migrations.CreateModel(
            name='AlignmentQC',
            fields=[
                ('dataprocess_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='osqpipe.DataProcess')),
                ('alignment', models.ForeignKey(to='osqpipe.Alignment', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'db_table': 'alignment_qc',
                'verbose_name': 'Alignment QC',
            },
            bases=('osqpipe.dataprocess',),
        ),
        migrations.CreateModel(
            name='AlnQCfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('filename', models.CharField(unique=True, max_length=1024)),
                ('checksum', models.CharField(max_length=128)),
                ('description', models.TextField(null=True, blank=True)),
                ('date', models.DateField(auto_now_add=True)),
                ('archive_date', models.DateField(null=True, blank=True)),
                ('alignmentqc', models.ForeignKey(to='osqpipe.AlignmentQC', on_delete=django.db.models.deletion.PROTECT)),
                ('archive', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, blank=True, to='osqpipe.ArchiveLocation', null=True)),
                ('filetype', models.ForeignKey(to='osqpipe.Filetype', on_delete=django.db.models.deletion.PROTECT)),
            ],
            options={
                'ordering': ['filename'],
                'db_table': 'alnqcfile',
                'verbose_name': 'Alignment QC file',
            },
        ),
    ]
