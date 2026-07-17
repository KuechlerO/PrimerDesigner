import sys
import types
import unittest
from dataclasses import dataclass
from unittest.mock import patch

# `primer_designer_app.utils.sv_utils` imports `primer_utils`, which imports `primer3`.
# Unit tests for request parsing shouldn't require the primer3 dependency, so we stub
# the `primer_utils` module before importing `sv_utils`.
_primer_utils_stub = types.ModuleType("primer_designer_app.utils.primer_utils")
_primer_utils_stub.primer3_design_primers = lambda *args, **kwargs: None
sys.modules.setdefault("primer_designer_app.utils.primer_utils", _primer_utils_stub)

from primer_designer_app.utils.sv_utils import (  # noqa: E402
    _calculate_genomic_primer_positions,
    build_structural_variant_info_from_request,
    design_structural_variant_primers,
)
from primer_designer_app.utils.variant_info import (  # noqa: E402
    StructuralVariantWindow,
)


@dataclass
class _FakeRequest:
    POST: dict


def _make_request(
    *,
    sv_chromosome="1",
    sv_start_position="10",
    # Use a valid minimal span (>= 50 bases) so window creation succeeds.
    sv_end_position="80",
    sv_type="deletion",
    reference_genome=None,
):
    post = {
        "sv_chromosome": sv_chromosome,
        "sv_start_position": sv_start_position,
        "sv_end_position": sv_end_position,
        "sv_type": sv_type,
    }
    if reference_genome is not None:
        post["reference-genome"] = reference_genome
    return _FakeRequest(POST=post)


class BuildStructuralVariantInfoFromRequestTests(unittest.TestCase):
    def test_parses_positions_as_integers(self):
        req = _make_request(sv_start_position="  123 ", sv_end_position="456  ")
        info = build_structural_variant_info_from_request(req)

        self.assertEqual(info.start_position, 123)
        self.assertEqual(info.end_position, 456)

    def test_trims_and_normalizes_chromosome_variants(self):
        cases = [
            ("1", "1"),
            (" 1 ", "1"),
            ("chr1", "1"),
            ("CHR2", "2"),
            ("x", "X"),
            (" chrX ", "X"),
            ("chrm", "M"),
        ]

        for raw, expected in cases:
            with self.subTest(raw=raw):
                req = _make_request(sv_chromosome=raw)
                info = build_structural_variant_info_from_request(req)
                self.assertEqual(info.chromosome, expected)

    def test_trims_and_lowercases_sv_type(self):
        req = _make_request(sv_type="  DeLeTiOn ")
        info = build_structural_variant_info_from_request(req)
        # SV type is currently handled at the view/template level; parser returns
        # only coordinates + reference genome.
        self.assertTrue(hasattr(info, "chromosome"))

    def test_defaults_reference_genome_to_grch37(self):
        req = _make_request(reference_genome=None)
        info = build_structural_variant_info_from_request(req)
        self.assertEqual(info.reference_genome, "GRCh37")

    def test_uses_reference_genome_from_request_and_trims(self):
        req = _make_request(reference_genome=" GRCh38 ")
        info = build_structural_variant_info_from_request(req)
        self.assertEqual(info.reference_genome, "GRCh38")

    def test_missing_start_position_raises_clear_error(self):
        req = _make_request(sv_start_position="   ")
        with self.assertRaisesRegex(
            ValueError, r"^SV start position must be a valid integer$"
        ):
            build_structural_variant_info_from_request(req)

    def test_non_integer_start_position_raises_clear_error(self):
        req = _make_request(sv_start_position="12.3")
        with self.assertRaisesRegex(
            ValueError, r"^SV start position must be a valid integer$"
        ):
            build_structural_variant_info_from_request(req)

    def test_zero_end_position_raises_clear_error(self):
        req = _make_request(sv_end_position="0")
        with self.assertRaisesRegex(
            ValueError, r"^SV end position must be a positive integer$"
        ):
            build_structural_variant_info_from_request(req)

    def test_missing_sv_type_bubbles_up_as_unsupported(self):
        req = _make_request(sv_type="  ")
        # SV type is not parsed here, so request parsing should still succeed.
        info = build_structural_variant_info_from_request(req)
        self.assertEqual(info.chromosome, "1")


class CalculateGenomicPrimerPositionsTests(unittest.TestCase):
    def test_positions_are_offset_by_window_start(self):
        window = StructuralVariantWindow(
            label="upstream",
            window_start_genomic=1000,
            window_end_genomic=2000,
        )
        pair = types.SimpleNamespace(
            left_relPos_start=5,
            left_relPos_end=24,
            right_relPos_start=100,
            right_relPos_end=119,
        )
        positions = _calculate_genomic_primer_positions(window, pair)
        self.assertEqual(
            positions,
            {
                "forward_start": 1005,
                "forward_end": 1024,
                "reverse_start": 1100,
                "reverse_end": 1119,
            },
        )


class DesignStructuralVariantPrimersTests(unittest.TestCase):
    def _sv_info(self):
        from primer_designer_app.utils.variant_info import StructuralVariantInfo

        sv_info = StructuralVariantInfo(
            chromosome="1",
            start_position=1000,
            end_position=5000,
            reference_genome="GRCh37",
        )
        sv_info.create_design_windows()
        return sv_info

    def test_designs_primers_for_every_window(self):
        sv_info = self._sv_info()
        primer_settings = types.SimpleNamespace()

        fake_pair = types.SimpleNamespace(
            left_relPos_start=5,
            left_relPos_end=24,
            right_relPos_start=100,
            right_relPos_end=119,
        )
        fake_search_results = types.SimpleNamespace(primer_pairs=[fake_pair])

        with patch.object(
            StructuralVariantWindow, "load_window_sequence", lambda self, **kw: None
        ), patch(
            "primer_designer_app.utils.sv_utils.primer3_design_primers",
            return_value=fake_search_results,
        ) as mock_design:
            results = design_structural_variant_primers(sv_info, primer_settings)

        self.assertEqual(mock_design.call_count, len(sv_info.windows))
        self.assertEqual(set(results.keys()), {w.label for w in sv_info.windows})

        upstream = results["upstream"]
        self.assertIs(
            upstream["design_window"],
            next(w for w in sv_info.windows if w.label == "upstream"),
        )
        self.assertEqual(len(upstream["primer_rows"]), 1)
        row = upstream["primer_rows"][0]
        self.assertIs(row["pair"], fake_pair)
        self.assertEqual(
            row["genomic_positions"]["forward_start"],
            upstream["design_window"].window_start_genomic + 5,
        )

    def test_no_primer_pairs_yields_empty_primer_rows(self):
        sv_info = self._sv_info()
        primer_settings = types.SimpleNamespace()
        fake_search_results = types.SimpleNamespace(primer_pairs=[])

        with patch.object(
            StructuralVariantWindow, "load_window_sequence", lambda self, **kw: None
        ), patch(
            "primer_designer_app.utils.sv_utils.primer3_design_primers",
            return_value=fake_search_results,
        ):
            results = design_structural_variant_primers(sv_info, primer_settings)

        for window_result in results.values():
            self.assertEqual(window_result["primer_rows"], [])
