from django.db import migrations, models


def _table_columns(schema_editor, table):
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("PRAGMA table_info(%s)" % schema_editor.quote_name(table))
        return {row[1] for row in cursor.fetchall()}


def add_snp_awareness_columns(apps, schema_editor):
    """Add SNP columns only if missing (schema may predate committed 0005)."""
    settings_table = "primer_designer_app_primersettingsmodel"
    summary_table = "primer_designer_app_designresultssummary"
    settings_cols = _table_columns(schema_editor, settings_table)
    summary_cols = _table_columns(schema_editor, summary_table)

    # #region agent log
    try:
        import json
        import time
        from pathlib import Path

        for log_path in (
            Path("/app/django_data/debug-migrate.ndjson"),
            Path(".cursor/debug-5ef680.log"),
        ):
            try:
                log_path.parent.mkdir(parents=True, exist_ok=True)
                with log_path.open("a") as fh:
                    fh.write(
                        json.dumps(
                            {
                                "sessionId": "5ef680",
                                "runId": "migrate-0005",
                                "hypothesisId": "H4",
                                "location": "0005_snp_awareness.py:add_snp_awareness_columns",
                                "message": "Idempotent 0005 column check",
                                "data": {
                                    "has_check_known_snps": "check_known_snps"
                                    in settings_cols,
                                    "has_snp_analysis_data": "snp_analysis_data"
                                    in summary_cols,
                                    "will_add_check_known_snps": "check_known_snps"
                                    not in settings_cols,
                                    "will_add_snp_analysis_data": "snp_analysis_data"
                                    not in summary_cols,
                                },
                                "timestamp": int(time.time() * 1000),
                            }
                        )
                        + "\n"
                    )
            except OSError:
                pass
    except Exception:
        pass
    # #endregion

    if "check_known_snps" not in settings_cols:
        schema_editor.execute(
            "ALTER TABLE %s ADD COLUMN check_known_snps bool NOT NULL DEFAULT 0"
            % schema_editor.quote_name(settings_table)
        )
    if "snp_analysis_data" not in summary_cols:
        schema_editor.execute(
            "ALTER TABLE %s ADD COLUMN snp_analysis_data text NULL"
            % schema_editor.quote_name(summary_table)
        )


class Migration(migrations.Migration):

    dependencies = [
        ("primer_designer_app", "0004_target_padding_remove_use_case"),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
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
            ],
            database_operations=[
                migrations.RunPython(
                    add_snp_awareness_columns, migrations.RunPython.noop
                ),
            ],
        ),
    ]
