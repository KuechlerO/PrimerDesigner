import json
import logging
import re
from typing import List, Optional, Sequence

from primer_designer_app.models import PrimerSettingsModel
from primer_designer_app.utils.variant_info import (
    AllelicVariantInfo,
    TranscriptVariantInfo,
    SequenceVariantInfo,
    GenomicVariantInfo,
    IndelType,
)
from primer_designer_app.utils.display_utils import (
    compute_display_bounds,
    shift_template_hits_for_display,
)
from primer_designer_app.utils.primer_utils import PrimerPairResult

LOGGER = logging.getLogger(__name__)


def vcf_hits_for_display(
    vcf_applied: List[dict], display_offset: int, display_length: int
) -> List[dict]:
    """Shift VCF template coordinates into the displayed sequence slice."""
    return shift_template_hits_for_display(
        list(vcf_applied or []), display_offset, display_length
    )


def vcf_hits_json_for_display(
    vcf_applied: List[dict], display_offset: int, display_length: int
) -> str:
    return json.dumps(vcf_hits_for_display(vcf_applied, display_offset, display_length))


def snp_hits_json_for_display(
    snp_hits: List[dict], display_offset: int, display_length: int
) -> str:
    return json.dumps(
        shift_template_hits_for_display(snp_hits, display_offset, display_length)
    )


def map_variant_content(var_info: AllelicVariantInfo) -> tuple[str, str]:
    """Maps variant info to indel type and notation for highlighting in HTML output."""
    if var_info.indel_type in [IndelType.NONE, IndelType.SNV]:
        LOGGER.debug(f"1")
        return ["", f"{var_info.ref_bases}>{var_info.new_bases}"]
    elif var_info.indel_type == IndelType.DEL:
        LOGGER.debug(f"2")
        return ["Del:", f"{var_info.ref_bases}/-"]
    elif var_info.indel_type == IndelType.INS:
        LOGGER.debug(f"3")
        return ["Ins:", f"-/{var_info.new_bases}"]
    elif var_info.indel_type == IndelType.DELINS:
        LOGGER.debug(f"4")
        return ["delins:", f"{var_info.ref_bases}/{var_info.new_bases}"]
    else:
        raise ValueError(
            f"Unsupported indel type for variant content mapping: {var_info.indel_type}"
        )


def create_hgvs_notation(var_info: AllelicVariantInfo) -> str:
    """
    Build HGVS-like notation for transcript or genomic variant info.
    """
    LOGGER.debug(f"Creating HGVS notation for variant: {var_info}")

    def get_final_hgvs_construct(
        var_info: AllelicVariantInfo, prefix, var_position
    ) -> str:
        LOGGER.debug(f"Variant position for HGVS notation: {var_position}")
        hgvs_notation = ""
        if var_info.indel_type in [IndelType.SNV]:
            var_pos_final = var_position[0]
            hgvs_notation = (
                f"{prefix}{var_pos_final}{var_info.ref_bases}>{var_info.new_bases}"
            )
        elif var_info.indel_type == IndelType.DEL:
            mut_start = var_position[0]
            mut_end = var_position[1]
            hgvs_notation = f"{prefix}{mut_start}_{mut_end}del"
        elif var_info.indel_type == IndelType.DELINS:
            mut_start = var_position[0]
            mut_end = var_position[1]
            hgvs_notation = f"{prefix}{mut_start}_{mut_end}delins{var_info.new_bases}"
        elif var_info.indel_type == IndelType.INS:
            mut_pos = var_position[0]
            hgvs_notation = f"{prefix}{mut_pos}_{mut_pos+1}ins{var_info.new_bases}"
        else:
            raise ValueError(
                f"Unsupported indel type for HGVS notation: {var_info.indel_type}"
            )
        return hgvs_notation + f" [{var_info.ref_genome}]"

    prefix = ""
    if isinstance(var_info, GenomicVariantInfo):
        var_position = var_info.get_genomic_pos()
        if var_info.gene_ID:
            prefix = f"{var_info.gene_symbol}({var_info.gene_ID}) "
        prefix_array = prefix + f"chr{var_info.genomic_pos.get('chr')}: g."
        return get_final_hgvs_construct(var_info, prefix_array, var_position)

    elif isinstance(var_info, TranscriptVariantInfo):
        var_info.get_genomic_pos()
        prefix_coding = f"{var_info.gene_symbol}({var_info.transcript_id}): c."
        prefix_genomic = f"chr{var_info.genomic_pos.get('chr')}: g."
        # return f"""{get_final_hgvs_construct(var_info, prefix_coding, relative_var_position)} - \
        #             {get_final_hgvs_construct(var_info, prefix_genomic, genomic_var_position)}"""
        return f"{get_final_hgvs_construct(var_info, prefix_coding, var_info.relative_pos)}"

    elif isinstance(var_info, SequenceVariantInfo):
        prefix = ""
        return get_final_hgvs_construct(var_info, prefix, var_info.relative_pos)
    else:
        raise ValueError(
            f"Unsupported AllelicVariantInfo type for HGVS notation: {type(var_info)}"
        )


def transform_rel_primer_pos(primerF_range, primerR_range, genomic_pos):
    for p_range in (primerF_range, primerR_range):
        p_range[0] += genomic_pos - 1000 - 1
        p_range[1] += genomic_pos - 1000 - 1
    return primerF_range, primerR_range


def _slice_variant_info_for_display(
    var_info: AllelicVariantInfo,
    target_start: int,
    target_len: int,
    primer_pairs: Sequence[PrimerPairResult],
) -> tuple[int, int, int, int, int, int, int, int]:
    """
    Slice var_info in place for UI display; return display offset and adjusted coordinates
    for the active primer pair (first in primer_pairs).
    """
    template_len = len(var_info.ref_seq)
    var_lo, var_hi = var_info.relative_pos
    pair = primer_pairs[0]
    primerF_start, primerF_end = pair.left_relPos_start, pair.left_relPos_end
    primerR_start, primerR_end = pair.right_relPos_start, pair.right_relPos_end

    display_start, display_end = compute_display_bounds(
        template_len, var_lo, var_hi, target_start, target_len, primer_pairs
    )

    saved_ref = var_info.ref_seq
    var_info.relative_pos
    saved_bases = var_info.ref_bases

    var_info.ref_seq = saved_ref[display_start:display_end]
    var_info.relative_pos = (var_lo - display_start, var_hi - display_start)
    if saved_bases:
        var_info.ref_bases = var_info.ref_seq[
            var_info.relative_pos[0] : var_info.relative_pos[1] + 1
        ]

    adj = display_start
    return (
        display_start,
        primerF_start - adj,
        primerF_end - adj,
        primerR_start - adj,
        primerR_end - adj,
        target_start - adj,
        display_end - display_start,
    )


def html_visualize_sequence(
    prim_settings: PrimerSettingsModel,
    var_info: AllelicVariantInfo,
    prim_pair: PrimerPairResult,
    *,
    all_primer_pairs: Optional[Sequence[PrimerPairResult]] = None,
) -> tuple[str, int, int, list[dict]]:
    """Return (HTML, display_offset, display_length, display_chunks) for the sequence view."""
    pairs = list(all_primer_pairs) if all_primer_pairs else [prim_pair]
    if prim_pair not in pairs:
        pairs = [prim_pair] + pairs

    saved_ref = var_info.ref_seq
    saved_rel = var_info.relative_pos
    saved_bases = var_info.ref_bases
    display_offset = 0
    display_length = len(var_info.ref_seq)
    highlighted_sequence = ""
    display_chunks: list[dict] = []
    orig_ref_bases = var_info.ref_bases or ""
    orig_new_bases = var_info.new_bases or ""

    try:
        (
            display_offset,
            primerF_start,
            primerF_end,
            primerR_start,
            primerR_end,
            target_start,
            display_length,
        ) = _slice_variant_info_for_display(
            var_info,
            prim_settings.target[0],
            prim_settings.target[1],
            pairs,
        )
        # Primer3 coordinates are on the mutated template; plain chunks must match.
        _slice_plain = var_info.get_seq("mutated")
        annotated_seq = var_info.get_seq("input")

        reg_expr = r"(\[[^\]]+\])"

        def _highlight_variant(match: re.Match) -> str:
            capture = match.group(1)
            var_text = capture.replace("[", "[" + var_info.indel_type.value + ":")
            return f"<span class='highlight-mutation'>{var_text}</span>"

        highlighted_sequence = re.sub(
            reg_expr, _highlight_variant, annotated_seq, count=1
        )

        display_chunks = build_allele_display_chunks(
            _slice_plain,
            var_info,
            allele="wt",
            width=100,
            ref_bases=orig_ref_bases,
            new_bases=orig_new_bases,
            highlight_snv_allele=False,
        )

        # Only highlight the variant itself (highlight-mutation).
        # The previous "highlight-target" wrapper around surrounding flank/target
        # area is intentionally removed to reduce visual noise.
    finally:
        var_info.ref_seq = saved_ref
        var_info.relative_pos = saved_rel
        var_info.ref_bases = saved_bases

    LOGGER.debug(
        "Display highlighted sequence length: %s (offset %s)",
        len(highlighted_sequence),
        display_offset,
    )
    return highlighted_sequence, display_offset, display_length, display_chunks


def _normalize_indel_type(var_info: AllelicVariantInfo) -> IndelType:
    t = getattr(var_info, "indel_type", None)
    if isinstance(t, IndelType):
        return t
    if isinstance(t, str):
        try:
            return IndelType(t)
        except ValueError:
            pass
    return IndelType.NONE


def _hgvs_input_on_plain(
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


def _allele_annotated_seq(
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
    indel_type = _normalize_indel_type(var_info)
    if allele == "mut" and indel_type == IndelType.DELINS:
        new_u = (new_bases or "").upper()
        if new_u:
            hi = lo + len(new_u) - 1
    return _hgvs_input_on_plain(
        plain, lo, hi, indel_type, ref_bases, new_bases, allele=allele
    )


def _highlight_annotated_variant(
    annotated_seq: str, var_info: AllelicVariantInfo
) -> str:
    """Wrap the first bracket group with highlight-mutation (same as SNV/Indel view)."""
    reg_expr = r"(\[[^\]]+\])"

    def _highlight_variant(match: re.Match) -> str:
        capture = match.group(1)
        indel_type = _normalize_indel_type(var_info)
        var_text = capture.replace("[", "[" + indel_type.value + ":")
        return f"<span class='highlight-mutation'>{var_text}</span>"

    return re.sub(reg_expr, _highlight_variant, annotated_seq, count=1)


def _highlight_allele_base_in_snv_bracket(
    html: str,
    var_info: AllelicVariantInfo,
    *,
    allele: str,
    ref_bases: str,
    new_bases: str,
) -> str:
    """
    For SNV, highlight the discriminating base inside the bracket:
    WT highlights ref base, MUT highlights alt base (notation stays [ref>alt]).
    """
    indel_type = _normalize_indel_type(var_info)
    if indel_type != IndelType.SNV:
        return html

    ref_bases_u = (ref_bases or "").upper()
    new_bases_u = (new_bases or "").upper()
    if not ref_bases_u or not new_bases_u:
        return html

    ref = ref_bases_u[0]
    alt = new_bases_u[0]
    chosen = ref if allele == "wt" else alt

    # Only touch the bracket that was already wrapped in highlight-mutation.
    pat = re.compile(
        r"(<span class='highlight-mutation'>\[[A-Z]+:)([ACGTN])>([ACGTN])(\]</span>)"
    )

    def _sub(m: re.Match) -> str:
        a = m.group(2)
        b = m.group(3)
        if chosen == a:
            a = f"<span class='highlight-primer'>{a}</span>"
        elif chosen == b:
            b = f"<span class='highlight-primer'>{b}</span>"
        return m.group(1) + a + ">" + b + m.group(4)

    return pat.sub(_sub, html, count=1)


def build_allele_display_chunks(
    plain_template: str,
    var_info: AllelicVariantInfo,
    *,
    allele: str,
    width: int,
    ref_bases: str,
    new_bases: str,
    highlight_snv_allele: bool = False,
) -> list[dict]:
    """
    Build paired plain/HTML chunks for AS-PCR display.

    Plain chunks align 1:1 with Primer3 template coordinates; HTML chunks add
    HGVS bracket notation and variant styling for the same template ranges.
    """
    if plain_template is None:
        return []
    try:
        width = int(width)
    except Exception:
        width = 50
    if width <= 0:
        width = 50

    lo, hi = var_info.relative_pos
    indel_type = _normalize_indel_type(var_info)
    chunks: list[dict] = []
    for i in range(0, len(plain_template), width):
        plain_chunk = plain_template[i : i + width]
        start_1 = i + 1
        chunk_end_0 = i + len(plain_chunk)
        if lo < chunk_end_0 and hi >= i:
            lo_local = lo - i
            hi_local = hi - i
            lo_local_used = lo_local
            hi_local_used = hi_local

            # Important: var_info.relative_pos spans the REF allele.
            # For DELINS, the MUT plain template has ALT length at this locus; if we remove
            # REF-length bases from the MUT template, downstream highlights shift by 1.
            if allele == "mut" and indel_type == IndelType.DELINS:
                new_u = (new_bases or "").upper()
                if new_u:
                    hi_local_used = lo_local_used + len(new_u) - 1
            annotated = _hgvs_input_on_plain(
                plain_chunk,
                lo_local_used,
                hi_local_used,
                indel_type,
                ref_bases,
                new_bases,
                allele=allele,
            )
            html_chunk = _highlight_annotated_variant(annotated, var_info)
            if highlight_snv_allele:
                html_chunk = _highlight_allele_base_in_snv_bracket(
                    html_chunk,
                    var_info,
                    allele=allele,
                    ref_bases=ref_bases,
                    new_bases=new_bases,
                )
        else:
            html_chunk = plain_chunk
        chunks.append({"start": start_1, "plain": plain_chunk, "html": html_chunk})
    return chunks


def html_visualize_sequence_allele_specific(
    prim_settings: PrimerSettingsModel,
    var_info: AllelicVariantInfo,
    prim_pair: PrimerPairResult,
    *,
    allele: str,
    all_primer_pairs: Optional[Sequence[PrimerPairResult]] = None,
) -> tuple[str, int, int, list[dict]]:
    """
    AS-PCR sequence view: HGVS bracket annotation on WT and MUT templates.

    Primer3 coordinates map to the plain template; JS highlightPrimerRegion skips
    letters inside ``[-/…]`` and SNV bracket groups (same rules for WT and MUT).
    """
    pairs = list(all_primer_pairs) if all_primer_pairs else [prim_pair]
    if prim_pair not in pairs:
        pairs = [prim_pair] + pairs

    saved_ref = var_info.ref_seq
    saved_rel = var_info.relative_pos
    saved_bases = var_info.ref_bases
    orig_ref_bases = var_info.ref_bases or ""
    orig_new_bases = var_info.new_bases or ""

    if allele == "mut":
        template_full = var_info.get_seq("mutated")
    else:
        template_full = var_info.ref_seq

    var_info.ref_seq = template_full
    display_offset = 0
    display_length = len(template_full)
    highlighted_sequence = ""
    display_chunks: list[dict] = []

    try:
        (
            display_offset,
            _primerF_start,
            _primerF_end,
            _primerR_start,
            _primerR_end,
            _target_start,
            display_length,
        ) = _slice_variant_info_for_display(
            var_info,
            prim_settings.target[0],
            prim_settings.target[1],
            pairs,
        )
        # var_info.ref_seq is already the sliced WT or MUT template (set above).
        # Do NOT call get_seq("mutated") here: that would re-apply DELINS/DEL using
        # ref-relative coordinates on an already-mutated slice and shift downstream
        # highlighting by (ref_len - alt_len).
        _slice_plain = var_info.ref_seq
        annotated = _allele_annotated_seq(
            var_info,
            ref_bases=orig_ref_bases,
            new_bases=orig_new_bases,
            allele=allele,
        )
        highlighted_sequence = _highlight_annotated_variant(annotated, var_info)
        highlighted_sequence = _highlight_allele_base_in_snv_bracket(
            highlighted_sequence,
            var_info,
            allele=allele,
            ref_bases=orig_ref_bases,
            new_bases=orig_new_bases,
        )
        display_chunks = build_allele_display_chunks(
            _slice_plain,
            var_info,
            allele=allele,
            width=50,
            ref_bases=orig_ref_bases,
            new_bases=orig_new_bases,
            highlight_snv_allele=True,
        )
    finally:
        var_info.ref_seq = saved_ref
        var_info.relative_pos = saved_rel
        var_info.ref_bases = saved_bases

    return highlighted_sequence, display_offset, display_length, display_chunks
