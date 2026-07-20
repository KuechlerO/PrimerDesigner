from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("primer_designer_app", "0002_primersettingsmodel_primer3_overrides"),
    ]

    operations = [
        migrations.AddField(
            model_name="primersettingsmodel",
            name="do_insilico_pcr",
            field=models.BooleanField(default=False),
        ),
    ]
