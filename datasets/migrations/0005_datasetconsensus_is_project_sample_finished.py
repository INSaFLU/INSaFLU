# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2022-06-07 23:43
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('datasets', '0004_consensus_number_of_locus'),
    ]

    operations = [
        migrations.AddField(
            model_name='datasetconsensus',
            name='is_project_sample_finished',
            field=models.BooleanField(default=False),
        ),
    ]
