# -*- coding: utf-8 -*-
# Generated by Django 1.11.29 on 2024-07-09 12:32
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pathogen_identification', '0050_auto_20240603_1533'),
    ]

    operations = [
        migrations.AddField(
            model_name='rawreferencecompoundmodel',
            name='screening_count',
            field=models.IntegerField(default=0),
        ),
    ]