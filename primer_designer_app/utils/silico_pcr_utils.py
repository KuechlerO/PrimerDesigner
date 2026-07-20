"""Helpers for the standalone SilicoPCR mode (user-supplied primers → Dicey)."""

from __future__ import annotations

import re
from types import SimpleNamespace

from primer_designer_app.exceptions import InvalidInputError
from primer_designer_app.utils.insilico_analysis import (
    do_insilico_analysis,
    insilico_reference_description,
)
from primer_designer_app.utils.primer_utils import PrimerPairResult

MIN_PRIMER_LEN = 10
MAX_PRIMER_LEN = 60
_PRIMER_ALPHABET_RE = re.compile(r"^[ACGTN]+$", re.IGNORECASE)


def normalize_primer_sequence(raw: str | None, field_name: str) -> str:
    """Strip whitespace and validate a user-supplied primer sequence."""
    if raw is None:
        raise InvalidInputError(f"{field_name} must not be empty.")
    seq = "".join(str(raw).split()).upper()
    if not seq:
        raise InvalidInputError(f"{field_name} must not be empty.")
    if not _PRIMER_ALPHABET_RE.match(seq):
        raise InvalidInputError(
            f"{field_name} must contain only DNA bases A, C, G, T, or N."
        )
    if len(seq) < MIN_PRIMER_LEN or len(seq) > MAX_PRIMER_LEN:
        raise InvalidInputError(
            f"{field_name} must be between {MIN_PRIMER_LEN} and "
            f"{MAX_PRIMER_LEN} bp (got {len(seq)} bp)."
        )
    return seq


def build_silico_settings(
    reference_genome: str, amplicon_check: str
) -> SimpleNamespace:
    """Minimal settings object for Dicey (no Primer3 fields required)."""
    genome = (reference_genome or "GRCh37").strip()
    if genome not in ("GRCh37", "GRCh38"):
        raise InvalidInputError(
            f"Unsupported reference genome: {genome}. Use GRCh37 or GRCh38."
        )
    check = (amplicon_check or "genome").strip().lower()
    if check not in ("genome", "transcriptome"):
        raise InvalidInputError(
            "Search context must be Genome or Transcriptome for SilicoPCR."
        )
    context = "genomic" if check == "genome" else "transcriptomic"
    return SimpleNamespace(
        reference_genome=genome,
        context=context,
        do_insilico_pcr=True,
    )


def build_primer_pair(forward: str, reverse: str) -> PrimerPairResult:
    """Wrap validated primer sequences as a PrimerPairResult for Dicey."""
    return PrimerPairResult(
        index=0,
        left_seq=forward,
        right_seq=reverse,
        penalty=0.0,
        product_size=0,
        left_relPos_start=None,
        left_relPos_end=None,
        right_relPos_start=None,
        right_relPos_end=None,
        amplicons=[],
        insilico_seq="",
        insilico_status=None,
        insilico_error_detail=None,
    )


def run_silico_pcr(
    *,
    forward_raw: str,
    reverse_raw: str,
    reference_genome: str,
    amplicon_check: str,
) -> tuple[PrimerPairResult, SimpleNamespace, str]:
    """
    Validate inputs, run Dicey, and return (pair, settings, reference_note).
    """
    forward = normalize_primer_sequence(forward_raw, "Forward primer")
    reverse = normalize_primer_sequence(reverse_raw, "Reverse primer")
    settings = build_silico_settings(reference_genome, amplicon_check)
    pair = build_primer_pair(forward, reverse)
    do_insilico_analysis(settings, [pair])
    note = insilico_reference_description(settings)
    return pair, settings, note
