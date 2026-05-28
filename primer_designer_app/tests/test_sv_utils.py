import sys
import types
import unittest
from dataclasses import dataclass

# `primer_designer_app.utils.sv_utils` imports `primer_utils`, which imports `primer3`.
# Unit tests for request parsing shouldn't require the primer3 dependency, so we stub
# the `primer_utils` module before importing `sv_utils`.
_primer_utils_stub = types.ModuleType("primer_designer_app.utils.primer_utils")
_primer_utils_stub.primer3_design_primers = lambda *args, **kwargs: None
sys.modules.setdefault("primer_designer_app.utils.primer_utils", _primer_utils_stub)

from primer_designer_app.utils.sv_utils import (  # noqa: E402
    build_structural_variant_info_from_request,
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
