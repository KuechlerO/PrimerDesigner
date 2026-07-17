"""Sequence display window helpers (no Django / primer3 dependency)."""

from typing import List, Protocol, Sequence, Union

DISPLAY_FLANK = 250
# Wider flank for DOCX reports so both primers and the full amplicon stay visible.
REPORT_DISPLAY_FLANK = 500
# Default characters per line when splitting annotated sequences for display.
DEFAULT_CHUNK_WIDTH = 100

GenomicCoord = Union[int, str, float, None]


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """True if inclusive integer intervals [a_start, a_end] and [b_start, b_end] overlap."""
    return a_start <= b_end and b_start <= a_end


def format_genomic_coord(value: GenomicCoord) -> str:
    """Format a genomic coordinate with period thousands separators (19.997.582)."""
    if value is None or value == "":
        return ""
    try:
        return f"{int(value):,}".replace(",", ".")
    except (TypeError, ValueError):
        return str(value)


def format_genomic_range(
    start: GenomicCoord,
    end: GenomicCoord,
    *,
    sep: str = " - ",
) -> str:
    """Format a start–end genomic range for display."""
    left = format_genomic_coord(start)
    right = format_genomic_coord(end)
    if not left and not right:
        return ""
    if not left:
        return right
    if not right:
        return left
    return f"{left}{sep}{right}"


class PrimerPairCoords(Protocol):
    left_relPos_start: int
    left_relPos_end: int
    right_relPos_start: int
    right_relPos_end: int


def _primer_pair_bounds(pair: PrimerPairCoords) -> tuple[int, int]:
    lo = min(
        pair.left_relPos_start,
        pair.left_relPos_end,
        pair.right_relPos_start,
        pair.right_relPos_end,
    )
    hi = max(
        pair.left_relPos_start,
        pair.left_relPos_end,
        pair.right_relPos_start,
        pair.right_relPos_end,
    )
    return lo, hi


def compute_display_bounds(
    template_len: int,
    var_lo: int,
    var_hi: int,
    target_start: int,
    target_len: int,
    primer_pairs: Sequence[PrimerPairCoords],
    *,
    flank: int = DISPLAY_FLANK,
) -> tuple[int, int]:
    """Return (display_start, display_end) half-open indices on the full template."""
    target_end = target_start + target_len - 1
    span_lo = min(var_lo, target_start)
    span_hi = max(var_hi, target_end)
    for pair in primer_pairs:
        p_lo, p_hi = _primer_pair_bounds(pair)
        span_lo = min(span_lo, p_lo)
        span_hi = max(span_hi, p_hi)

    display_start = max(0, span_lo - flank)
    display_end = min(template_len, span_hi + flank + 1)
    return display_start, display_end


def compute_report_display_bounds(
    template_len: int,
    var_lo: int,
    var_hi: int,
    target_start: int,
    target_len: int,
    primer_pair: PrimerPairCoords,
) -> tuple[int, int]:
    """Report sequence window: same span as the UI but with a wider flank."""
    return compute_display_bounds(
        template_len,
        var_lo,
        var_hi,
        target_start,
        target_len,
        [primer_pair],
        flank=REPORT_DISPLAY_FLANK,
    )


def shift_template_hits_for_display(
    hits: List[dict], display_offset: int, display_length: int
) -> List[dict]:
    shifted: List[dict] = []
    for row in hits or []:
        ts = int(row["template_start"]) - display_offset
        te = int(row["template_end"]) - display_offset
        if te < 0 or ts >= display_length:
            continue
        shifted.append(
            {
                **row,
                "template_start": max(0, ts),
                "template_end": min(display_length - 1, te),
            }
        )
    return shifted
