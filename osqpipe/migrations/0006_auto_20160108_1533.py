# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models

from django.contrib.auth.models import User
from rest_framework.authtoken.models import Token

'''
Creates initial token set for all users.
'''

def create_user_tokens(apps, schema_editor):

    for user in User.objects.all():
        Token.objects.get_or_create(user=user)

class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0005_auto_20151218_1551'),
        ('authtoken', '0001_initial'),
    ]

    operations = [
        migrations.RunPython(create_user_tokens),
    ]
