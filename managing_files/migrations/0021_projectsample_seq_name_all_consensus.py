# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2021-04-17 14:17
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('managing_files', '0020_auto_20210414_2101'),
    ]

    operations = [
        migrations.AddField(
            model_name='projectsample',
            name='seq_name_all_consensus',
            field=models.CharField(blank=True, max_length=200, null=True),
        ),
    ]
