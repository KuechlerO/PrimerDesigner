from primer_designer_app.utils.variant_info import StructuralVariantInfo
from primer_designer_app.utils.primer_utils import primer3_design_primers


def _parse_positive_integer(raw_value: str, field_name: str) -> int:
    raw_value = str(raw_value).strip()
    if raw_value == '':
        raise ValueError(f"{field_name} is required")

    try:
        value = int(raw_value)
    except ValueError as exc:
        raise ValueError(f"{field_name} must be an integer") from exc

    if value < 1:
        raise ValueError(f"{field_name} must be greater than 0")

    return value


def build_structural_variant_info_from_request(request) -> StructuralVariantInfo:
    chromosome = request.POST.get('sv_chromosome', '').strip()
    start_position = _parse_positive_integer(
        request.POST.get('sv_start_position', ''),
        'SV start position',
    )
    end_position = _parse_positive_integer(
        request.POST.get('sv_end_position', ''),
        'SV end position',
    )
    structural_variant_type = request.POST.get('sv_type', '').strip()
    reference_genome = request.POST.get('reference-genome', 'GRCh37').strip()

    return StructuralVariantInfo(
        chromosome=chromosome,
        start_position=start_position,
        end_position=end_position,
        structural_variant_type=structural_variant_type,
        reference_genome=reference_genome,
    )


def calculate_sv_genomic_primer_positions(primer_pair, design_window: StructuralVariantInfo) -> dict:
    if primer_pair is None:
        return {
            'forward_start': None,
            'forward_end': None,
            'reverse_start': None,
            'reverse_end': None,
        }

    if primer_pair.left_relPos_start is None or primer_pair.left_relPos_end is None:
        forward_start = None
        forward_end = None
    else:
        forward_start = design_window.window_start_genomic + primer_pair.left_relPos_start
        forward_end = design_window.window_start_genomic + primer_pair.left_relPos_end

    if primer_pair.right_relPos_start is None or primer_pair.right_relPos_end is None:
        reverse_start = None
        reverse_end = None
    else:
        reverse_start = design_window.window_start_genomic + primer_pair.right_relPos_start
        reverse_end = design_window.window_start_genomic + primer_pair.right_relPos_end

    return {
        'forward_start': forward_start,
        'forward_end': forward_end,
        'reverse_start': reverse_start,
        'reverse_end': reverse_end,
    }


def design_structural_variant_primers(
    structural_variant_info: StructuralVariantInfo,
    primer_settings,
) -> dict:
    design_windows = structural_variant_info.create_design_windows()

    sv_results = {}

    for design_window in design_windows:
        design_window.prepare_for_primer_design()

        primer_settings.set_target(design_window.get_target_interval_in_window())

        primer_search_results = primer3_design_primers(
            primer_settings,
            design_window,
        )

        best_primer_pair = None
        if primer_search_results.primer_pairs:
            best_primer_pair = primer_search_results.primer_pairs[0]

        genomic_positions = calculate_sv_genomic_primer_positions(
            best_primer_pair,
            design_window,
        )

        sv_results[design_window.label] = {
            'design_window': design_window,
            'best_primer_pair': best_primer_pair,
            'genomic_positions': genomic_positions,
            'all_primer_pairs': primer_search_results.primer_pairs,
        }

    return sv_results
