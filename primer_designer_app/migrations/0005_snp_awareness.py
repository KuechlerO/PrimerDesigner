from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("primer_designer_app", "0004_target_padding_remove_use_case"),
    ]

    operations = [
        migrations.AddField(
            model_name="primersettingsmodel",
            name="check_known_snps",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="designresultssummary",
            name="snp_analysis_data",
            field=models.JSONField(blank=True, null=True),
        ),
    ]
