# -*- coding: utf-8 -*-
# Generated by Django 1.11.29 on 2024-02-05 12:27
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pathogen_identification', '0028_teleflusample'),
    ]

    operations = [
        migrations.AlterField(
            model_name='metareference',
            name='description',
            field=models.CharField(blank=True, max_length=200, null=True),
        ),
    ]
