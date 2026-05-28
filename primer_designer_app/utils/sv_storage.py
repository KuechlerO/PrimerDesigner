"""Serialize / deserialize structural variant design results (no Django imports)."""

from primer_designer_app.utils.variant_info import StructuralVariantInfo

SV_DESIGN_TYPE = "structural_variant"
SV_WINDOW_ORDER = ("upstream", "internal_1", "internal_2", "downstream")


def serialize_structural_variant_info(
    structural_variant_info: StructuralVariantInfo,
) -> dict:
    """JSON-safe SV query metadata (no window sequences)."""
    return {
        "design_type": SV_DESIGN_TYPE,
        "chromosome": structural_variant_info.chromosome,
        "start_position": structural_variant_info.start_position,
        "end_position": structural_variant_info.end_position,
        "reference_genome": structural_variant_info.reference_genome,
        "windows": [
            {
                "label": window.label,
                "window_start_genomic": window.window_start_genomic,
                "window_end_genomic": window.window_end_genomic,
            }
            for window in structural_variant_info.windows
        ],
    }


def serialize_sv_results_for_storage(sv_results: dict) -> dict:
    """Convert in-memory SV primer results to JSON-storable dict."""
    stored = {}
    for label, result in sv_results.items():
        design_window = result["design_window"]
        stored[label] = {
            "window": {
                "label": design_window.label,
                "window_start_genomic": design_window.window_start_genomic,
                "window_end_genomic": design_window.window_end_genomic,
            },
            "primer_rows": [
                {
                    "pair": row["pair"].to_dict(),
                    "genomic_positions": row["genomic_positions"],
                }
                for row in result["primer_rows"]
            ],
        }
    return stored


def deserialize_sv_results_from_storage(stored: dict, primer_pair_from_dict) -> dict:
    """Rebuild SV results dict with PrimerPairResult objects."""
    windows = stored.get("windows") or {}
    results = {}
    for label in SV_WINDOW_ORDER:
        if label not in windows:
            continue
        window_data = windows[label]
        primer_rows = [
            {
                "pair": primer_pair_from_dict(row["pair"]),
                "genomic_positions": row["genomic_positions"],
            }
            for row in window_data.get("primer_rows", [])
        ]
        results[label] = {
            "window": window_data["window"],
            "primer_rows": primer_rows,
        }
    for label, window_data in windows.items():
        if label in results:
            continue
        primer_rows = [
            {
                "pair": primer_pair_from_dict(row["pair"]),
                "genomic_positions": row["genomic_positions"],
            }
            for row in window_data.get("primer_rows", [])
        ]
        results[label] = {
            "window": window_data["window"],
            "primer_rows": primer_rows,
        }
    return results
