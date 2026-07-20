from django.db import migrations, models


def use_case_to_target_padding(apps, schema_editor):
    PrimerSettings = apps.get_model("primer_designer_app", "PrimerSettingsModel")
    for row in PrimerSettings.objects.all():
        if getattr(row, "use_case", None) == "qPCR":
            row.target_padding = 30
        else:
            row.target_padding = 50
        row.save(update_fields=["target_padding"])


def target_padding_to_use_case(apps, schema_editor):
    PrimerSettings = apps.get_model("primer_designer_app", "PrimerSettingsModel")
    for row in PrimerSettings.objects.all():
        pad = getattr(row, "target_padding", 50) or 50
        row.use_case = "qPCR" if pad == 30 else "PCR"
        row.save(update_fields=["use_case"])


class Migration(migrations.Migration):

    dependencies = [
        ("primer_designer_app", "0003_primersettingsmodel_do_insilico_pcr"),
    ]

    operations = [
        migrations.AddField(
            model_name="primersettingsmodel",
            name="target_padding",
            field=models.IntegerField(default=50),
        ),
        migrations.RunPython(use_case_to_target_padding, target_padding_to_use_case),
        migrations.RemoveField(
            model_name="primersettingsmodel",
            name="use_case",
        ),
    ]
