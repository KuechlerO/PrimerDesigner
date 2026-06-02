"""Validation helpers for primer design inputs and Primer3 outcomes."""

from __future__ import annotations

from typing import Any

from primer_designer_app.exceptions import (
    InvalidReferenceSequenceError,
    NoPrimerPairsFoundError,
)

# Ensembl soft-masking (mask_feature=1) replaces intron/UTR bases with N.
_MASKED_FRACTION_THRESHOLD = 0.95


def masked_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    upper = seq.upper()
    return sum(1 for base in upper if base == "N") / len(upper)


def validate_reference_sequence_for_design(seq: str) -> None:
    """Reject templates that cannot support primer design (empty or almost all N)."""
    if not seq or not seq.strip():
        raise InvalidReferenceSequenceError(
            "No reference sequence was retrieved for the requested region."
        )

    if masked_fraction(seq) >= _MASKED_FRACTION_THRESHOLD:
        raise InvalidReferenceSequenceError(
            "The reference window is almost entirely masked (N bases). "
            "If this is not expected, please check the reference sequence and "
            "make sure that chromosome and position match the selected genome build "
            "(GRCh37 vs GRCh38)."
        )


def validate_primer_search_results(
    results: Any,
    *,
    context_label: str = "design",
) -> None:
    if not getattr(results, "primer_pairs", None):
        raise NoPrimerPairsFoundError(
            f"Primer3 did not return any primer pairs for this {context_label}. "
            "Try relaxing primer length, Tm, or GC limits, widening the product "
            "size range, or adjusting target padding."
        )
