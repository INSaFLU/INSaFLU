# -*- coding: utf-8 -*-
# Generated by Django 1.11.29 on 2023-11-21 18:33
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pathogen_identification', '0017_auto_20230905_1405'),
    ]

    operations = [
        migrations.AddField(
            model_name='finalreport',
            name='error_rate',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='finalreport',
            name='quality_avg',
            field=models.FloatField(blank=True, null=True),
        ),
    ]