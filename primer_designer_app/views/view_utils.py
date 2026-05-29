# View utilities for primer design app
import logging
from enum import Enum
from typing import Tuple


from primer_designer_app.models import PrimerSettingsModel, DesignResultsSummary
from primer_designer_app.utils.variant_info import (
    AllelicVariantInfo,
    TranscriptVariantInfo,
    GenomicVariantInfo,
    SequenceVariantInfo,
    ReferenceType,
    VARIANT_FLANKING,
)
from primer_designer_app.utils.primer3_post import (
    P3_FORM_DEFAULTS,
    PRIMER3_OVERRIDE_FIELDS,
    parse_primer3_overrides_from_post,
)
from primer_designer_app.utils.primer_utils import primer3_design_primers
from primer_designer_app.utils.vcf_utils import parse_vcf_upload

from primer_designer_app.exceptions import (
    InvalidTranscriptVersionError,
    InvalidTranscriptIdError,
    InvalidTranscriptInputError,
    ExonExonJunctionError,
)

LOGGER = logging.getLogger(__name__)


def _process_genome_pos_snv_input(input_pos: str, end_offset=0) -> dict:
    """Parse human input like 'chr1:12345' into genomic_pos dict.

    Args:
        input_pos (str): User input string representing genome position.
        end_offset (int): Offset to add to the end position.

    Returns:
        dict: Parsed genomic position with keys 'chr', 'pos', and 'strand_type'.
    """

    chromosome = ""
    position = [-1, -1]

    cleaned = input_pos.replace(",", "").strip()
    # simple parser: accept chrN:POS or chrMT:POS
    if ":" not in cleaned:
        raise ValueError(
            f"Invalid genomic position '{input_pos}'. "
            "Expected format ChrN:position (e.g. chr13:2655000)."
        )
    left, right = cleaned.split(":", 1)
    chromosome = left.lower().replace("chr", "").upper()
    try:
        position = int(right)
    except ValueError:
        raise ValueError(f"Invalid position in input: {position}")

    return {
        "chr": chromosome,
        "pos": [position, position + end_offset],
        "strand_type": "sense",
    }


def _process_genome_pos_indel_input(
    input_chr: str, input_start: str, input_end: str
) -> dict:
    """Parse human input like 'chr1:12345-12350' into genomic_pos dict.
    Args:
        input_chr (str):    User input string representing chromosome.
        input_start (str):  User input string representing start position.
        input_end (str):    User input string representing end position.
    Returns:
        dict: Parsed genomic position with keys 'chr', 'pos', and 'strand_type'.
    """
    chromosome = input_chr.lower().replace("chr", "").upper()
    try:
        start_pos = int(input_start)
        end_pos = int(input_end)
    except ValueError:
        raise ValueError(f"Invalid position in input: {input_start} or {input_end}")

    assert (
        start_pos <= end_pos
    ), "Start position must be less than or equal to end position."
    return {"chr": chromosome, "pos": [start_pos, end_pos], "strand_type": "sense"}


def _get_post(request, name, default=""):
    return request.POST.get(name, default)


def _parse_snp_check(request) -> bool:
    return _get_post(request, "snp-check", "").lower() in ("on", "1", "true", "yes")


def build_form_data_from_request(request, **extra) -> dict:
    """
    Build a template context dict so primer settings and optional page fields
    repopulate after POST (avoids resetting dialog inputs to HTML defaults).
    """
    data = {
        "reference_genome": _get_post(request, "reference-genome", "GRCh37"),
        "amplicon_check": _get_post(request, "amplicon-check", "none"),
        "snp_check": "on" if _parse_snp_check(request) else "",
        "tm": _get_post(request, "tm", "60"),
        "gc_content": _get_post(request, "gc_content", "50"),
        "product_size_min": _get_post(request, "product_size_min", "400"),
        "product_size_max": _get_post(request, "product_size_max", "800"),
        "target_padding": _get_post(request, "target_padding", "50"),
        "max_poly_X": _get_post(request, "max_poly_X", "4"),
    }
    for key, _kind in PRIMER3_OVERRIDE_FIELDS:
        field_name = f"p3_{key}"
        data[field_name] = _get_post(
            request, field_name, P3_FORM_DEFAULTS.get(field_name, "")
        )
    data.update(extra)
    return data


def _parse_amplicon_check(request) -> Tuple[bool, str]:
    """
    Single UI control: none | genome | transcriptome.
    Maps to do_insilico_pcr and PrimerSettings.context for Dicey.
    """
    v = _get_post(request, "amplicon-check", "none")
    if v not in ("none", "genome", "transcriptome"):
        v = "none"
    if v == "none":
        return False, "genomic"
    if v == "genome":
        return True, "genomic"
    return True, "transcriptomic"


def _parse_optional_vcf_upload(request, chromosome: str):
    """Return parsed VCF records for the design chromosome, or None if no upload."""
    upload = request.FILES.get("vcf_file")
    if not upload or not getattr(upload, "name", ""):
        return None
    try:
        return parse_vcf_upload(upload, chromosome)
    except ValueError as exc:
        raise ValueError(f"VCF upload: {exc}") from exc


def _parse_target_padding(request) -> int:
    try:
        v = int(_get_post(request, "target_padding", "50"))
    except ValueError:
        return 50
    return max(1, min(500, v))


def build_primer_settings(request) -> PrimerSettingsModel:
    do_insilico, context = _parse_amplicon_check(request)
    return PrimerSettingsModel(
        target_padding=_parse_target_padding(request),
        tm=int(_get_post(request, "tm", "60")),
        gc=int(_get_post(request, "gc_content", "50")),
        reference_genome=_get_post(request, "reference-genome", "GRCh37"),
        productsize_range=[
            int(_get_post(request, "product_size_min", "400")),
            int(_get_post(request, "product_size_max", "800")),
        ],
        max_poly_x=int(_get_post(request, "max_poly_X", "4")),
        primer3_overrides=parse_primer3_overrides_from_post(request),
        do_insilico_pcr=do_insilico,
        context=context,
        check_known_snps=_parse_snp_check(request),
    )


def _build_variant_info(request, input_type: str) -> AllelicVariantInfo:
    ref_genome = _get_post(request, "reference-genome", "GRCh37")

    if input_type == "genomic_snv":
        genomic_pos = _get_post(request, "genom_pos", "")
        new_base = _get_post(request, "new_base", "").strip().upper()
        if len(new_base) != 1 or new_base not in "ACGT":
            raise ValueError(
                "For genomic SNV input, provide exactly one alternate base (A, C, G, or T) "
                "in the New base field."
            )
        gpos = _process_genome_pos_snv_input(genomic_pos, len(new_base) - 1)
        vcf_records = _parse_optional_vcf_upload(request, gpos["chr"])
        rel_pos = (
            None
            if vcf_records
            else (VARIANT_FLANKING, VARIANT_FLANKING + len(new_base) - 1)
        )
        variant_info = GenomicVariantInfo(
            genomic_pos=gpos,
            new_bases=new_base,
            ref_genome=ref_genome,
            relative_pos=rel_pos or (0, 0),
            vcf_records=vcf_records,
        )

    elif input_type == "genomic_indel":
        indelChrom = _get_post(request, "IndelChrom", "")
        indelStart = _get_post(request, "IndelStart", "")
        indelEnd = _get_post(request, "IndelEnd", "")
        indelIns = _get_post(request, "IndelIns", "")
        if indelIns.isnumeric():
            indelIns = "N" * int(indelIns)  # Convert numeric input to string of Ns
        gpos = _process_genome_pos_indel_input(indelChrom, indelStart, indelEnd)
        vcf_records = _parse_optional_vcf_upload(request, gpos["chr"])
        indel_span = int(indelEnd) - int(indelStart)
        rel_pos = (
            None
            if vcf_records
            else (VARIANT_FLANKING, VARIANT_FLANKING + indel_span - 1)
        )
        variant_info = GenomicVariantInfo(
            genomic_pos=gpos,
            new_bases=indelIns,
            ref_genome=ref_genome,
            relative_pos=rel_pos or (0, 0),
            vcf_records=vcf_records,
        )

    elif input_type in ["transcript_snv", "transcript_indel"]:
        transcript_id = _get_post(request, "Transcript-ID", "")
        post_ref_type = _get_post(request, "Reference", "")
        if post_ref_type == "cdna":
            reference_type = ReferenceType.CDNA
        elif post_ref_type == "cds":
            reference_type = ReferenceType.CDS
        else:
            raise ValueError(f"Unknown reference type: {post_ref_type}")

        if input_type == "transcript_snv":
            position = _get_post(request, "Position", "")
            new_bases = _get_post(request, "IDnew_base", "")
            relative_pos = [int(position) - 1, int(position) - 1]
            LOGGER.debug(
                f"Parsed transcript SNV input: position={position}, new_bases={new_bases}, relative_pos={relative_pos}"
            )

        elif input_type == "transcript_indel":
            indel_start = _get_post(request, "IdIndelStart", "")
            indel_end = _get_post(request, "IdIndelEnd", "")
            new_bases = _get_post(request, "IdIndelIns", "")
            if new_bases.isnumeric():
                new_bases = "N" * int(
                    new_bases
                )  # Convert numeric input to string of Ns
            relative_pos = [int(indel_start) - 1, int(indel_end) - 1]

        variant_info = TranscriptVariantInfo(
            ref_genome=ref_genome,
            transcript_id=transcript_id,
            reference_type=reference_type,
            new_bases=new_bases,
            relative_pos=relative_pos,
        )

    elif input_type == "sequence_input":
        input_seq = _get_post(request, "sequence", "")
        variant_info = SequenceVariantInfo(
            input_seq=input_seq,
            ref_genome=ref_genome,
        )
    else:
        raise ValueError(f"Unknown input_type: {input_type}")

    LOGGER.debug(f"Built AllelicVariantInfo: {variant_info}")
    return variant_info


def handle_genomic_snv(request, primer_settings: PrimerSettingsModel):
    variantInfo = _build_variant_info(request, "genomic_snv")
    primer_settings.set_target(variantInfo.relative_pos)
    return _design_primers_and_return_searchID(variantInfo, primer_settings)


def handle_genomic_indel(request, primer_settings):
    # Replace INDEL input with HGVS and VCF input format
    variantInfo = _build_variant_info(request, "genomic_indel")
    primer_settings.set_target(variantInfo.relative_pos)
    return _design_primers_and_return_searchID(variantInfo, primer_settings)


# TODO: Check for invalid transcript input & invalid letter input
def handle_transcript_input(request, primer_settings):
    # --- 1. Create AllelicVariantInfo from transcript input ---
    if _get_post(request, "Position", ""):
        variantInfo = _build_variant_info(request, "transcript_snv")
    elif _get_post(request, "IdIndelStart", "") and _get_post(
        request, "IdIndelEnd", ""
    ):
        variantInfo = _build_variant_info(request, "transcript_indel")
    else:
        # TODO: handle error properly
        raise InvalidTranscriptInputError(
            "The transcript input is incomplete or invalid."
        )
    # Check input validity
    # TODO: check
    # type int means that the Indel is only within and therefore only has one mapping/ one range
    if type(variantInfo.genomic_pos["pos"][0]) != int:
        raise ExonExonJunctionError("The variant affects an exon-exon junction.")

    # --- 2. Design primers (context / insilico flags come from build_primer_settings) ---
    primer_settings.set_target(variantInfo.relative_pos)
    LOGGER.debug(
        "Transcript input: context=%s do_insilico_pcr=%s target=%s",
        primer_settings.context,
        primer_settings.do_insilico_pcr,
        primer_settings.target,
    )

    return _design_primers_and_return_searchID(variantInfo, primer_settings)


def handle_sequence_input(request, primer_settings):
    variantInfo = _build_variant_info(request, "sequence_input")
    primer_settings.set_target(variantInfo.relative_pos)
    return _design_primers_and_return_searchID(variantInfo, primer_settings)


def handle_allele_specific_input(request, primer_settings: PrimerSettingsModel):
    """
    Allele-specific PCR (AS-PCR) entrypoint.

    Builds the same AllelicVariantInfo objects as SNV/Indel mode, but stores
    allele-specific Primer3 output (WT + MUT) in primer_search_results.
    """
    from primer_designer_app.utils.primer_utils import primer3_design_allele_specific

    # Determine input type using the same rules as SNV/Indel index
    if _get_post(request, "Transcript-ID", None):
        if _get_post(request, "Position", ""):
            variantInfo = _build_variant_info(request, "transcript_snv")
        elif _get_post(request, "IdIndelStart", "") and _get_post(
            request, "IdIndelEnd", ""
        ):
            variantInfo = _build_variant_info(request, "transcript_indel")
        else:
            raise InvalidTranscriptInputError(
                "The transcript input is incomplete or invalid."
            )
    elif _get_post(request, "genom_pos", None):
        variantInfo = _build_variant_info(request, "genomic_snv")
    elif _get_post(request, "IndelChrom", None):
        variantInfo = _build_variant_info(request, "genomic_indel")
    else:
        variantInfo = _build_variant_info(request, "sequence_input")

    primer_settings.set_target(variantInfo.relative_pos)
    primer_settings.do_insilico_pcr = False
    primer_settings.check_known_snps = False

    allele_specific_data = primer3_design_allele_specific(primer_settings, variantInfo)

    # Persist results in DesignResultsSummary (custom serialization)
    result_sum_obj = DesignResultsSummary()
    result_sum_obj.primer_settings = primer_settings

    # Serialize AllelicVariantInfo (copy pattern from save_primer_results)
    variant_info_dict = variantInfo.__dict__.copy()
    if isinstance(variant_info_dict.get("indel_type"), Enum):
        variant_info_dict["indel_type"] = variant_info_dict["indel_type"].value
    if isinstance(variant_info_dict.get("reference_type"), Enum):
        variant_info_dict["reference_type"] = variant_info_dict["reference_type"].value
    variant_info_dict["design_type"] = "allele_specific"
    result_sum_obj.variant_info_data = variant_info_dict

    result_sum_obj.primer_search_results = allele_specific_data
    result_sum_obj.snp_analysis_data = {}
    result_sum_obj.save()
    return result_sum_obj.id


def _design_primers_and_return_searchID(variant_info, primer_settings):
    """Final step: Design primers and generate output."""
    from primer_designer_app.utils.snp_awareness import (
        annotate_primer_pairs_with_snp_awareness,
    )

    LOGGER.debug(f"Context for primer design: {primer_settings.context}")
    primer_res = primer3_design_primers(primer_settings, variant_info)
    LOGGER.debug(f"Context 2 for primer design: {primer_settings.context}")

    target = primer_settings.target
    primer_target = (int(target[0]), int(target[1])) if target else None
    snp_analysis = annotate_primer_pairs_with_snp_awareness(
        variant_info,
        primer_res.primer_pairs,
        primer_settings.reference_genome,
        enabled=getattr(primer_settings, "check_known_snps", False),
        primer_target=primer_target,
    )

    result_sum_obj = DesignResultsSummary()
    result_sum_obj.save_primer_results(
        primer_res, primer_settings, variant_info, snp_analysis_data=snp_analysis
    )
    return result_sum_obj.id
