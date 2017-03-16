# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0010_alignmentqc_alnqcfile'),
    ]

    operations = [
        migrations.AddField(
            model_name='library',
            name='external_records',
            field=models.ManyToManyField(related_name='libraries', db_table=b'library_external_record', to='osqpipe.ExternalRecord'),
        ),
        migrations.AddField(
            model_name='sample',
            name='external_records',
            field=models.ManyToManyField(related_name='samples', db_table=b'sample_external_record', to='osqpipe.ExternalRecord'),
        ),
    ]
