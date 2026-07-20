#!/usr/bin/env python
"""Pre-migrate probe for Docker deploy debugging (session 5ef680)."""

from __future__ import annotations

import json
import os
import sqlite3
import sys
import time
from pathlib import Path


def _write_log(entry: dict) -> None:
    paths = [
        Path("/app/django_data/debug-migrate.ndjson"),
        Path(".cursor/debug-5ef680.log"),
    ]
    line = json.dumps(entry) + "\n"
    for path in paths:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("a") as fh:
                fh.write(line)
        except OSError:
            pass
    print(line, end="", file=sys.stderr)


def main() -> int:
    db_path = Path(os.environ.get("DJANGO_DB_PATH", "/app/django_data/db.sqlite3"))
    # Also resolve relative project path used in settings
    alt = Path("django_data/db.sqlite3")
    if not db_path.exists() and alt.exists():
        db_path = alt

    data = {
        "db_path": str(db_path),
        "db_exists": db_path.exists(),
        "db_size": db_path.stat().st_size if db_path.exists() else None,
        "has_snp_analysis_data": None,
        "has_check_known_snps": None,
        "applied_migrations": [],
        "has_0005_recorded": None,
    }

    # #region agent log
    _write_log(
        {
            "sessionId": "5ef680",
            "runId": os.environ.get("DEBUG_RUN_ID", "pre-migrate"),
            "hypothesisId": "H1",
            "location": "scripts/debug_migrate_probe.py",
            "message": "Runtime DB presence before migrate",
            "data": {
                "db_path": data["db_path"],
                "db_exists": data["db_exists"],
                "db_size": data["db_size"],
            },
            "timestamp": int(time.time() * 1000),
        }
    )
    # #endregion

    if not db_path.exists():
        _write_log(
            {
                "sessionId": "5ef680",
                "runId": os.environ.get("DEBUG_RUN_ID", "pre-migrate"),
                "hypothesisId": "H5",
                "location": "scripts/debug_migrate_probe.py",
                "message": "No DB yet — fresh migrate expected",
                "data": data,
                "timestamp": int(time.time() * 1000),
            }
        )
        return 0

    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute("PRAGMA table_info(primer_designer_app_designresultssummary)")
    summary_cols = {row[1] for row in cur.fetchall()}
    cur.execute("PRAGMA table_info(primer_designer_app_primersettingsmodel)")
    settings_cols = {row[1] for row in cur.fetchall()}
    cur.execute(
        "SELECT name FROM django_migrations WHERE app='primer_designer_app' ORDER BY id"
    )
    migs = [row[0] for row in cur.fetchall()]
    conn.close()

    data["has_snp_analysis_data"] = "snp_analysis_data" in summary_cols
    data["has_check_known_snps"] = "check_known_snps" in settings_cols
    data["applied_migrations"] = migs
    data["has_0005_recorded"] = any(m.startswith("0005") for m in migs)
    data["duplicate_column_risk"] = bool(
        data["has_snp_analysis_data"] and not data["has_0005_recorded"]
    )

    # #region agent log
    _write_log(
        {
            "sessionId": "5ef680",
            "runId": os.environ.get("DEBUG_RUN_ID", "pre-migrate"),
            "hypothesisId": "H2",
            "location": "scripts/debug_migrate_probe.py",
            "message": "Schema vs migration history",
            "data": data,
            "timestamp": int(time.time() * 1000),
        }
    )
    # #endregion

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
