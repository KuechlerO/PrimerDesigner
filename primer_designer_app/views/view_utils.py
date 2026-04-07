# View utilities for primer design app
import logging
from typing import Tuple


from primer_designer_app.models import PrimerSettingsModel, DesignResultsSummary
from primer_designer_app.utils.variant_info import (
    VariantInfo,
    TranscriptVariantInfo,
    GenomicVariantInfo,
    SequenceVariantInfo,
    ReferenceType,
    VARIANT_FLANKING,
)
from primer_designer_app.utils.primer3_post import parse_primer3_overrides_from_post
from primer_designer_app.utils.primer_utils import primer3_design_primers

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

    chromosome = ''
    position = [-1, -1]

    cleaned = input_pos.replace(',', '').strip()
    # simple parser: accept chrN:POS or chrMT:POS
    if ':' not in cleaned:
        LOGGER.debug('Invalid genome position input: %s', input_pos)
        return {'chr': chromosome, 'pos': position, 'strand_type': 'sense'}
    left, right = cleaned.split(':', 1)
    chromosome = left.lower().replace('chr', '').upper()
    try:
        position = int(right)
    except ValueError:
        raise ValueError(f"Invalid position in input: {position}")

    return {
        'chr': chromosome,
        'pos': [position, position + end_offset],
        'strand_type': 'sense',
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
    chromosome = input_chr.lower().replace('chr', '').upper()
    try:
        start_pos = int(input_start)
        end_pos = int(input_end)
    except ValueError:
        raise ValueError(f"Invalid position in input: {input_start} or {input_end}")

    assert (
        start_pos <= end_pos
    ), 'Start position must be less than or equal to end position.'
    return {'chr': chromosome, 'pos': [start_pos, end_pos], 'strand_type': 'sense'}


def _get_post(request, name, default=''):
    return request.POST.get(name, default)


def _parse_amplicon_check(request) -> Tuple[bool, str]:
    """
    Single UI control: none | genome | transcriptome.
    Maps to do_insilico_pcr and PrimerSettings.context for Dicey.
    """
    v = _get_post(request, 'amplicon-check', 'none')
    if v not in ('none', 'genome', 'transcriptome'):
        v = 'none'
    if v == 'none':
        return False, 'genomic'
    if v == 'genome':
        return True, 'genomic'
    return True, 'transcriptomic'


def build_primer_settings(request) -> PrimerSettingsModel:
    do_insilico, context = _parse_amplicon_check(request)
    return PrimerSettingsModel(
        use_case=_get_post(request, 'usecase', ''),
        tm=int(_get_post(request, 'tm', '60')),
        gc=int(_get_post(request, 'gc_content', '50')),
        reference_genome=_get_post(request, 'reference-genome', 'GRCh37'),
        productsize_range=[
            int(_get_post(request, 'product_size_min', '400')),
            int(_get_post(request, 'product_size_max', '800')),
        ],
        max_poly_x=int(_get_post(request, 'max_poly_X', '4')),
        primer3_overrides=parse_primer3_overrides_from_post(request),
        do_insilico_pcr=do_insilico,
        context=context,
    )


def _build_variant_info(request, input_type: str) -> VariantInfo:
    ref_genome = _get_post(request, 'reference-genome', 'GRCh37')

    if input_type == 'genomic_snv':
        genomic_pos = _get_post(request, 'genom_pos', '')
        new_base = _get_post(request, 'new_base', '')
        assert len(new_base) == 1, 'For SNV input, new_base must be a single character.'
        variant_info = GenomicVariantInfo(
            genomic_pos=_process_genome_pos_snv_input(genomic_pos, len(new_base) - 1),
            new_bases=new_base,
            ref_genome=ref_genome,
            relative_pos=(VARIANT_FLANKING, VARIANT_FLANKING + len(new_base) - 1),
        )

    elif input_type == 'genomic_indel':
        indelChrom = _get_post(request, 'IndelChrom', '')
        indelStart = _get_post(request, 'IndelStart', '')
        indelEnd = _get_post(request, 'IndelEnd', '')
        indelIns = _get_post(request, 'IndelIns', '')
        if indelIns.isnumeric():
            indelIns = 'N' * int(indelIns)  # Convert numeric input to string of Ns
        variant_info = GenomicVariantInfo(
            genomic_pos=_process_genome_pos_indel_input(
                indelChrom, indelStart, indelEnd
            ),
            new_bases=indelIns,
            ref_genome=ref_genome,
            relative_pos=(
                VARIANT_FLANKING,
                VARIANT_FLANKING + (int(indelEnd) - int(indelStart)) - 1,
            ),
        )

    elif input_type in ['transcript_snv', 'transcript_indel']:
        transcript_id = _get_post(request, 'Transcript-ID', '')
        post_ref_type = _get_post(request, 'Reference', '')
        if post_ref_type == 'cdna':
            reference_type = ReferenceType.CDNA
        elif post_ref_type == 'cds':
            reference_type = ReferenceType.CDS
        else:
            raise ValueError(f"Unknown reference type: {post_ref_type}")

        if input_type == 'transcript_snv':
            position = _get_post(request, 'Position', '')
            new_bases = _get_post(request, 'IDnew_base', '')
            relative_pos = [int(position) - 1, int(position) - 1]
            LOGGER.debug(
                f"Parsed transcript SNV input: position={position}, new_bases={new_bases}, relative_pos={relative_pos}"
            )

        elif input_type == 'transcript_indel':
            indel_start = _get_post(request, 'IdIndelStart', '')
            indel_end = _get_post(request, 'IdIndelEnd', '')
            new_bases = _get_post(request, 'IdIndelIns', '')
            if new_bases.isnumeric():
                new_bases = 'N' * int(
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

    elif input_type == 'sequence_input':
        input_seq = _get_post(request, 'sequence', '')
        variant_info = SequenceVariantInfo(
            input_seq=input_seq,
            ref_genome=ref_genome,
        )
    else:
        raise ValueError(f"Unknown input_type: {input_type}")

    LOGGER.debug(f"Built VariantInfo: {variant_info}")
    return variant_info


def handle_genomic_snv(request, primer_settings: PrimerSettingsModel):
    variantInfo = _build_variant_info(request, 'genomic_snv')
    primer_settings.set_target(variantInfo.relative_pos)
    return _design_primers_and_return_searchID(variantInfo, primer_settings)


def handle_genomic_indel(request, primer_settings):
    # Replace INDEL input with HGVS and VCF input format
    variantInfo = _build_variant_info(request, 'genomic_indel')
    primer_settings.set_target(variantInfo.relative_pos)
    return _design_primers_and_return_searchID(variantInfo, primer_settings)


# TODO: Check for invalid transcript input & invalid letter input
def handle_transcript_input(request, primer_settings):
    # --- 1. Create VariantInfo from transcript input ---
    if _get_post(request, 'Position', ''):
        variantInfo = _build_variant_info(request, 'transcript_snv')
    elif _get_post(request, 'IdIndelStart', '') and _get_post(
        request, 'IdIndelEnd', ''
    ):
        variantInfo = _build_variant_info(request, 'transcript_indel')
    else:
        # TODO: handle error properly
        raise InvalidTranscriptInputError(
            'The transcript input is incomplete or invalid.'
        )
    # Check input validity
    # TODO: check
    # type int means that the Indel is only within and therefore only has one mapping/ one range
    if type(variantInfo.genomic_pos['pos'][0]) != int:
        raise ExonExonJunctionError('The variant affects an exon-exon junction.')

    # --- 2. Design primers (context / insilico flags come from build_primer_settings) ---
    primer_settings.set_target(variantInfo.relative_pos)
    LOGGER.debug(
        'Transcript input: context=%s do_insilico_pcr=%s target=%s',
        primer_settings.context,
        primer_settings.do_insilico_pcr,
        primer_settings.target,
    )

    return _design_primers_and_return_searchID(variantInfo, primer_settings)


def handle_sequence_input(request, primer_settings):
    variantInfo = _build_variant_info(request, 'sequence_input')
    primer_settings.set_target(variantInfo.relative_pos)
    return _design_primers_and_return_searchID(variantInfo, primer_settings)


def _design_primers_and_return_searchID(variant_info, primer_settings):
    """Final step: Design primers and generate output."""
    LOGGER.debug(f"Context for primer design: {primer_settings.context}")
    primer_res = primer3_design_primers(primer_settings, variant_info)
    LOGGER.debug(f"Context 2 for primer design: {primer_settings.context}")

    result_sum_obj = DesignResultsSummary()
    result_sum_obj.save_primer_results(primer_res, primer_settings, variant_info)
    return result_sum_obj.id
