"""
Parse optional Primer3 global_args from POST (advanced Primer3 keys).

Field names in HTML: p3_<PRIMER_KEY>. Non-empty values are merged after base args in primer_utils.
"""

import logging
from typing import Any

LOGGER = logging.getLogger(__name__)

# (primer3_key, python_type_name)
# Defaults for repopulating the primer-parameters dialog after POST (HTML name -> value).
P3_FORM_DEFAULTS = {
    "p3_PRIMER_OPT_SIZE": "20",
    "p3_PRIMER_MIN_SIZE": "18",
    "p3_PRIMER_MAX_SIZE": "22",
    "p3_PRIMER_MIN_TM": "",
    "p3_PRIMER_MAX_TM": "",
    "p3_PRIMER_MIN_GC": "20",
    "p3_PRIMER_MAX_GC": "80",
    "p3_PRIMER_GC_CLAMP": "1",
    "p3_PRIMER_SALT_MONOVALENT": "50",
    "p3_PRIMER_DNA_CONC": "50",
    "p3_PRIMER_MAX_NS_ACCEPTED": "0",
    "p3_PRIMER_MAX_SELF_ANY": "12",
    "p3_PRIMER_MAX_SELF_END": "8",
    "p3_PRIMER_PAIR_MAX_COMPL_ANY": "12",
    "p3_PRIMER_PAIR_MAX_COMPL_END": "8",
    "p3_PRIMER_INSIDE_PENALTY": "1",
    "p3_PRIMER_INTERNAL_MAX_SELF_END": "8",
    "p3_PRIMER_INTERNAL_MAX_POLY_X": "100",
}

PRIMER3_OVERRIDE_FIELDS = [
    ("PRIMER_OPT_SIZE", "int"),
    ("PRIMER_MIN_SIZE", "int"),
    ("PRIMER_MAX_SIZE", "int"),
    ("PRIMER_MIN_TM", "float"),
    ("PRIMER_MAX_TM", "float"),
    ("PRIMER_MIN_GC", "float"),
    ("PRIMER_MAX_GC", "float"),
    ("PRIMER_GC_CLAMP", "int"),
    ("PRIMER_SALT_MONOVALENT", "float"),
    ("PRIMER_DNA_CONC", "float"),
    ("PRIMER_MAX_NS_ACCEPTED", "int"),
    ("PRIMER_MAX_SELF_ANY", "int"),
    ("PRIMER_MAX_SELF_END", "int"),
    ("PRIMER_PAIR_MAX_COMPL_ANY", "int"),
    ("PRIMER_PAIR_MAX_COMPL_END", "int"),
    ("PRIMER_INSIDE_PENALTY", "float"),
    ("PRIMER_INTERNAL_MAX_SELF_END", "int"),
    ("PRIMER_INTERNAL_MAX_POLY_X", "int"),
]


def _coerce(raw: str, kind: str) -> Any:
    raw = raw.strip()
    if kind == "int":
        return int(raw)
    if kind == "float":
        return float(raw)
    raise ValueError(kind)


def parse_primer3_overrides_from_post(request) -> dict:
    """Return a dict of Primer3 keys for non-empty p3_* POST fields (advanced overrides)."""
    out: dict[str, Any] = {}
    for key, kind in PRIMER3_OVERRIDE_FIELDS:
        raw = request.POST.get(f"p3_{key}")
        if raw is None or str(raw).strip() == "":
            continue
        try:
            out[key] = _coerce(str(raw), kind)
        except ValueError as e:
            LOGGER.warning("Skip invalid primer3 POST %s=%r: %s", key, raw, e)
    return out
