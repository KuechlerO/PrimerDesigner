import logging
import re

from primer_designer_app.models import PrimerSettingsModel
from primer_designer_app.utils.variant_info import (
    TranscriptVariantInfo,
    SequenceVariantInfo,
    VariantInfo,
    GenomicVariantInfo,
    IndelType,
)
from primer_designer_app.utils.primer_utils import PrimerPairResult


LOGGER = logging.getLogger(__name__)


def map_variant_content(var_info: VariantInfo) -> tuple[str, str]:
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


def create_hgvs_notation(var_info: VariantInfo) -> str:
    """
    Build HGVS-like notation for transcript or genomic variant info.
    """
    LOGGER.debug(f"Creating HGVS notation for variant: {var_info}")

    def get_final_hgvs_construct(var_info: VariantInfo, prefix, var_position) -> str:
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
        genomic_var_position = var_info.get_genomic_pos()
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
            f"Unsupported VariantInfo type for HGVS notation: {type(var_info)}"
        )


def transform_rel_primer_pos(primerF_range, primerR_range, genomic_pos):
    for p_range in (primerF_range, primerR_range):
        p_range[0] += genomic_pos - 1000 - 1
        p_range[1] += genomic_pos - 1000 - 1
    return primerF_range, primerR_range


def html_visualize_sequence(
    prim_settings: PrimerSettingsModel,
    var_info: VariantInfo,
    prim_pair: PrimerPairResult,
) -> str:
    """Generates an HTML string with the target region highlighted, including the variant and primer binding sites."""
    annotated_seq = var_info.get_seq("input")

    primerF_start, primerF_end = [
        prim_pair.left_relPos_start,
        prim_pair.left_relPos_end,
    ]
    primerR_start, primerR_end = [
        prim_pair.right_relPos_start,
        prim_pair.right_relPos_end,
    ]

    # 1. Get variant text -> extract text surrounded by [] for highlighting
    reg_expr = r"(\[[^\]]+\])"
    capture = re.findall(reg_expr, annotated_seq)[0]
    # Insert IndelType annotation
    var_text = capture.replace("[", "[" + var_info.indel_type.value + ":")

    start_target_region = max(prim_settings.target[0], primerF_end + 1)
    end_target_region = min(
        prim_settings.target[0] + prim_settings.target[1], primerR_start - 1
    )
    seq_target_region = annotated_seq[
        start_target_region:end_target_region
    ]  # TODO: Really?
    upstream_ref_pos = var_info.relative_pos[0] - start_target_region
    # Position downstream of variant in reference sequence
    downstream_ref_pos = upstream_ref_pos + len(capture)

    # 2. Create highlighted target region with highlighted variant
    seq_target_highlight = (
        f"<span class='highlight-target'>"  # Highlighting of target region
        + seq_target_region[:upstream_ref_pos]
        + f"<span class='highlight-mutation'>{var_text}</span>"  # Notation of variant
        + seq_target_region[downstream_ref_pos:]
        + "</span>"
    )

    highlighted_sequence = (
        annotated_seq[:start_target_region]
        + seq_target_highlight
        + annotated_seq[end_target_region:]
    )

    LOGGER.debug(f"Full highlighted sequence: {highlighted_sequence}")
    return highlighted_sequence
