from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pathogen_identification', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='referencesource',
            name='lineage_path',
            field=models.CharField(blank=True, max_length=1000, null=True),
        ),
    ]
