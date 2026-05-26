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
) -> tuple[str, int, int]:
    """Return (HTML, display_offset, display_length) for the sequence view."""
    pairs = list(all_primer_pairs) if all_primer_pairs else [prim_pair]
    if prim_pair not in pairs:
        pairs = [prim_pair] + pairs

    saved_ref = var_info.ref_seq
    saved_rel = var_info.relative_pos
    saved_bases = var_info.ref_bases
    display_offset = 0
    display_length = len(var_info.ref_seq)

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
        annotated_seq = var_info.get_seq("input")

        reg_expr = r"(\[[^\]]+\])"

        def _highlight_variant(match: re.Match) -> str:
            capture = match.group(1)
            var_text = capture.replace("[", "[" + var_info.indel_type.value + ":")
            return f"<span class='highlight-mutation'>{var_text}</span>"

        highlighted_sequence = re.sub(
            reg_expr, _highlight_variant, annotated_seq, count=1
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
    return highlighted_sequence, display_offset, display_length
