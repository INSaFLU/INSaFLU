# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2021-09-19 16:27
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('managing_files', '0022_project_number_passed_sequences'),
    ]

    operations = [
        migrations.AddField(
            model_name='sample',
            name='date_deleted_processed_fastq',
            field=models.DateTimeField(blank=True, null=True, verbose_name='Date attached'),
        ),
        migrations.AddField(
            model_name='sample',
            name='is_deleted_processed_fastq',
            field=models.BooleanField(default=False),
        ),
    ]
