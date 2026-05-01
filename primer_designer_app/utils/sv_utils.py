from copy import deepcopy

from primer_designer_app.utils.variant_info import (
    StructuralVariantInfo,
    StructuralVariantWindow,
)
from primer_designer_app.utils.primer_utils import primer3_design_primers


def _normalize_chromosome(chromosome: str) -> str:
    return chromosome.lower().replace("chr", "").upper()


def _parse_positive_integer(value: str, field_name: str) -> int:
    try:
        parsed_value = int(str(value).strip())
    except (TypeError, ValueError):
        raise ValueError(f"{field_name} must be a valid integer")

    if parsed_value <= 0:
        raise ValueError(f"{field_name} must be a positive integer")

    return parsed_value


def build_structural_variant_info_from_request(request) -> StructuralVariantInfo:
    chromosome = _normalize_chromosome(request.POST.get("sv_chromosome", "").strip())
    start_position = _parse_positive_integer(
        request.POST.get("sv_start_position", ""),
        "SV start position",
    )
    end_position = _parse_positive_integer(
        request.POST.get("sv_end_position", ""),
        "SV end position",
    )
    reference_genome = request.POST.get("reference-genome", "GRCh37").strip()

    if not chromosome:
        raise ValueError("Chromosome must not be empty")

    structural_variant_info = StructuralVariantInfo(
        chromosome=chromosome,
        start_position=start_position,
        end_position=end_position,
        reference_genome=reference_genome,
    )

    structural_variant_info.create_design_windows()
    return structural_variant_info


def _calculate_genomic_primer_positions(
    design_window: StructuralVariantWindow,
    primer_pair,
) -> dict:
    return {
        "forward_start": design_window.window_start_genomic
        + primer_pair.left_relPos_start,
        "forward_end": design_window.window_start_genomic + primer_pair.left_relPos_end,
        "reverse_start": design_window.window_start_genomic
        + primer_pair.right_relPos_start,
        "reverse_end": design_window.window_start_genomic
        + primer_pair.right_relPos_end,
    }


def design_structural_variant_primers(
    structural_variant_info: StructuralVariantInfo,
    primer_settings,
) -> dict:
    results = {}

    for design_window in structural_variant_info.windows:
        design_window.load_window_sequence(
            chromosome=structural_variant_info.chromosome,
            reference_genome=structural_variant_info.reference_genome,
        )

        design_window.set_default_target(target_length=150)

        window_primer_settings = deepcopy(primer_settings)
        window_primer_settings.target = design_window.get_primer3_target()

        primer_search_results = primer3_design_primers(
            window_primer_settings,
            design_window,
        )

        primer_rows = []
        for pair in primer_search_results.primer_pairs:
            genomic_positions = _calculate_genomic_primer_positions(
                design_window,
                pair,
            )
            primer_rows.append(
                {
                    "pair": pair,
                    "genomic_positions": genomic_positions,
                }
            )

        results[design_window.label] = {
            "design_window": design_window,
            "primer_rows": primer_rows,
        }

    return results
