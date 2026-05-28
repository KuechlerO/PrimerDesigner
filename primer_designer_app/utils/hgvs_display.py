"""
UI-oriented HGVS-style bracket annotation helpers.

These functions build bracketed variant notation (e.g. [G>A], [-/ATGC], [ATGC/-],
[GCAGG/GGTC]) on top of a plain template sequence. They are used by the web UI and
DOCX report rendering.
"""

from __future__ import annotations

from primer_designer_app.utils.variant_info import AllelicVariantInfo, IndelType


def normalize_indel_type(var_info: AllelicVariantInfo) -> IndelType:
    """Return a normalized IndelType from variant objects."""

    t = getattr(var_info, "indel_type", None)
    if isinstance(t, IndelType):
        return t
    if isinstance(t, str):
        try:
            return IndelType(t)
        except ValueError:
            pass
    return IndelType.NONE


def hgvs_input_on_plain(
    plain: str,
    lo: int,
    hi: int,
    indel_type: IndelType,
    ref_bases: str,
    new_bases: str,
    *,
    allele: str = "wt",
) -> str:
    """Apply HGVS input bracket notation on a plain template (WT or MUT)."""

    new_u = (new_bases or "").upper()
    alt_len = len(new_u) if new_u else 0
    # On MUT templates, inserted/replaced bases are already present in the plain sequence.
    # If we append `plain[lo:]` after the bracket, the ALT can appear twice.
    mut_suffix_start = lo + alt_len if allele == "mut" and alt_len else lo

    if indel_type == IndelType.SNV:
        return plain[:lo] + "[" + ref_bases + ">" + new_bases + "]" + plain[lo + 1 :]
    if indel_type == IndelType.INS:
        return plain[:lo] + "[-/" + new_bases + "]" + plain[mut_suffix_start:]
    if indel_type == IndelType.DEL:
        return plain[:lo] + "[" + ref_bases + "/-]" + plain[hi + 1 :]
    if indel_type == IndelType.DELINS:
        return (
            plain[:lo]
            + "["
            + ref_bases
            + "/"
            + new_bases
            + "]"
            + plain[mut_suffix_start if allele == "mut" and alt_len else hi + 1 :]
        )
    return plain


def allele_annotated_seq(
    var_info: AllelicVariantInfo,
    *,
    ref_bases: str,
    new_bases: str,
    allele: str = "wt",
) -> str:
    """
    Same HGVS input bracket notation on WT and MUT display templates.

    Uses original ref/alt alleles from variant metadata (not bases read back from
    a sliced mutated template, which would turn SNV [G>A] into [A>A]).
    """

    plain = var_info.ref_seq
    lo, hi = var_info.relative_pos
    indel_type = normalize_indel_type(var_info)
    if allele == "mut" and indel_type == IndelType.DELINS:
        new_u = (new_bases or "").upper()
        if new_u:
            hi = lo + len(new_u) - 1
    return hgvs_input_on_plain(
        plain, lo, hi, indel_type, ref_bases, new_bases, allele=allele
    )


def template_bases_consumed_by_bracket(
    body: str, plain_from_here: str, *, allele: str
) -> int:
    """
    Return how many template bases the bracket *consumes* at this locus.

    Mirrors the web UI behavior used for coordinate mapping around bracket notation.
    """

    b = body
    colon = b.find(":")
    if colon >= 0:
        b = b[colon + 1 :]
    b = b.strip()

    if ">" in b:
        return 1

    # INS: WT shows inserted bases which are not on the template. MUT may include
    # inserted bases in the template (sequence input); if so, consume them.
    if b.startswith("-/"):
        ins = b[2:]
        if ins and plain_from_here.upper().startswith(ins.upper()):
            return len(ins)
        return 0

    # DEL: WT consumes REF bases, MUT consumes 0.
    if b.endswith("/-"):
        ref = b[: b.index("/-")]
        return 0 if allele == "mut" else len(ref)

    # DELINS: prefer whichever side matches the current template prefix.
    if "/" in b:
        ref, alt = b.split("/", 1)
        if alt and plain_from_here.upper().startswith(alt.upper()):
            return len(alt)
        if ref and plain_from_here.upper().startswith(ref.upper()):
            return len(ref)
        return len(ref)

    return 1
