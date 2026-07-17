from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("primer_designer_app", "0001_initial"),
    ]

    operations = [
        migrations.AddField(
            model_name="primersettingsmodel",
            name="primer3_overrides",
            field=models.JSONField(blank=True, default=dict),
        ),
    ]
