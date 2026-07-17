import sys
import types
import unittest
from unittest.mock import patch

_primer_utils_stub = types.ModuleType("primer_designer_app.utils.primer_utils")
sys.modules.setdefault("primer_designer_app.utils.primer_utils", _primer_utils_stub)

from dataclasses import dataclass
from typing import Optional


@dataclass
class PrimerPairResult:
    index: int
    left_seq: str
    right_seq: str
    penalty: float
    product_size: int
    left_relPos_start: Optional[int] = None
    left_relPos_end: Optional[int] = None
    right_relPos_start: Optional[int] = None
    right_relPos_end: Optional[int] = None
    snp_status: Optional[str] = None
    snp_conflicts: Optional[list] = None


_primer_utils_stub.PrimerPairResult = PrimerPairResult
from primer_designer_app.utils.snp_awareness import (
    COMMON_VARIANT_MAF_THRESHOLD,
    SNP_STATUS_CAUTION,
    SNP_STATUS_CONFLICT,
    SNP_STATUS_ERROR,
    SNP_STATUS_NONE,
    SNP_STATUS_SKIPPED,
    _attach_maf_and_filter_common_variants,
    _classify_pair,
    _intervals_overlap,
    annotate_primer_pairs_with_snp_awareness,
    effective_maf,
    get_design_region_genomic,
)
from primer_designer_app.utils.variant_info import SequenceVariantInfo
from types import SimpleNamespace


class SnpAwarenessTests(unittest.TestCase):
    def test_snp_query_region_not_full_vcf_fetch_window(self):
        """SNP overlap must not use the entire VCF-expanded ref_seq span."""
        from types import SimpleNamespace

        var_info = SimpleNamespace(
            genomic_pos={"chr": "5", "pos": [80083380]},
            sequence_region_start=80034457,
            ref_seq="A" * 137238,
            relative_pos=(48923, 48923),
        )
        pair = PrimerPairResult(
            index=0,
            left_seq="A" * 20,
            right_seq="T" * 20,
            penalty=0.1,
            product_size=400,
            left_relPos_start=49002,
            left_relPos_end=49021,
            right_relPos_start=48987,
            right_relPos_end=49006,
        )
        region = get_design_region_genomic(
            var_info,
            primer_pairs=[pair],
            primer_target=(48923 - 50, 101),
        )
        self.assertIsNotNone(region)
        query_span = region["query_end"] - region["query_start"] + 1
        template_span = region["region_end"] - region["region_start"] + 1
        display_span = region["display_end"] - region["display_start"]
        self.assertLess(query_span, 2500)
        self.assertGreater(template_span, 100_000)
        self.assertGreater(region["query_start"], region["region_start"])
        self.assertLess(region["query_end"], region["region_end"])
        self.assertEqual(query_span, display_span)

    def test_effective_maf_from_top_level(self):
        self.assertAlmostEqual(
            effective_maf({"MAF": 0.174254}),
            0.174254,
        )

    def test_effective_maf_from_gnomad_global(self):
        variation = {
            "minor_allele": "G",
            "populations": [
                {"population": "gnomADg:ALL", "allele": "T", "frequency": 0.83},
                {"population": "gnomADg:ALL", "allele": "G", "frequency": 0.17},
            ],
        }
        self.assertAlmostEqual(effective_maf(variation), 0.17)

    def test_effective_maf_ignores_sample_level_pops(self):
        variation = {
            "populations": [
                {"population": "ALFA:SAMN1", "allele": "A", "frequency": 1.0},
            ],
        }
        self.assertIsNone(effective_maf(variation))

    def test_filter_common_variants_by_maf(self):
        hits = [
            {"id": "rs_common", "template_start": 1, "template_end": 1},
            {"id": "rs_rare", "template_start": 2, "template_end": 2},
        ]
        details = {
            "rs_common": {"MAF": 0.05},
            "rs_rare": {"MAF": 0.005},
        }
        filtered = _attach_maf_and_filter_common_variants(
            hits,
            details,
            maf_threshold=COMMON_VARIANT_MAF_THRESHOLD,
        )
        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered[0]["id"], "rs_common")
        self.assertEqual(filtered[0]["maf"], 0.05)

    def test_intervals_overlap(self):
        self.assertTrue(_intervals_overlap(5, 10, 8, 12))
        self.assertFalse(_intervals_overlap(5, 10, 11, 15))

    def test_classify_conflict_on_forward_binding(self):
        pair = PrimerPairResult(
            index=0,
            left_seq="A" * 20,
            right_seq="T" * 20,
            penalty=0.1,
            product_size=100,
            left_relPos_start=10,
            left_relPos_end=29,
            right_relPos_start=80,
            right_relPos_end=99,
        )
        hits = [
            {
                "id": "rs1",
                "template_start": 15,
                "template_end": 15,
                "genomic_start": 1,
                "genomic_end": 1,
                "alleles": "A/G",
            }
        ]
        status, details = _classify_pair(pair, hits)
        self.assertEqual(status, SNP_STATUS_CONFLICT)
        self.assertEqual(len(details), 1)

    def test_classify_caution_in_amplicon_only(self):
        pair = PrimerPairResult(
            index=0,
            left_seq="A" * 20,
            right_seq="T" * 20,
            penalty=0.1,
            product_size=100,
            left_relPos_start=10,
            left_relPos_end=19,
            right_relPos_start=50,
            right_relPos_end=69,
        )
        hits = [
            {
                "id": "rs2",
                "template_start": 25,
                "template_end": 25,
                "genomic_start": 1,
                "genomic_end": 1,
                "alleles": "C/T",
            }
        ]
        status, _ = _classify_pair(pair, hits)
        self.assertEqual(status, SNP_STATUS_CAUTION)

    def test_classify_none_when_no_overlap(self):
        pair = PrimerPairResult(
            index=0,
            left_seq="A" * 20,
            right_seq="T" * 20,
            penalty=0.1,
            product_size=100,
            left_relPos_start=10,
            left_relPos_end=19,
            right_relPos_start=80,
            right_relPos_end=99,
        )
        hits = [
            {
                "id": "rs3",
                "template_start": 200,
                "template_end": 200,
                "genomic_start": 1,
                "genomic_end": 1,
                "alleles": "G/A",
            }
        ]
        status, details = _classify_pair(pair, hits)
        self.assertEqual(status, SNP_STATUS_NONE)
        self.assertEqual(details, [])


def _genomic_var_info():
    return SimpleNamespace(
        genomic_pos={"chr": "7", "pos": [2000]},
        sequence_region_start=1000,
        ref_seq="A" * 3000,
        relative_pos=(999, 999),
    )


def _genomic_pair():
    return PrimerPairResult(
        index=0,
        left_seq="A" * 20,
        right_seq="T" * 20,
        penalty=0.1,
        product_size=200,
        left_relPos_start=900,
        left_relPos_end=919,
        right_relPos_start=1080,
        right_relPos_end=1099,
    )


class AnnotatePrimerPairsWithSnpAwarenessTests(unittest.TestCase):
    def test_disabled_marks_pairs_skipped(self):
        pair = _genomic_pair()
        summary = annotate_primer_pairs_with_snp_awareness(
            _genomic_var_info(), [pair], "GRCh37", enabled=False
        )
        self.assertEqual(summary["status"], SNP_STATUS_SKIPPED)
        self.assertEqual(summary["message"], "SNP check was not requested.")
        self.assertEqual(pair.snp_status, SNP_STATUS_SKIPPED)
        self.assertEqual(pair.snp_conflicts, [])

    def test_sequence_only_input_is_skipped(self):
        seq_var = SequenceVariantInfo(input_seq="AC[G>A]GT", ref_genome="GRCh37")
        self.assertIsNone(getattr(seq_var, "genomic_pos", None))
        pair = _genomic_pair()
        summary = annotate_primer_pairs_with_snp_awareness(
            seq_var, [pair], "GRCh37", enabled=True
        )
        self.assertEqual(summary["status"], SNP_STATUS_SKIPPED)
        self.assertIn("requires genomic or transcript coordinates", summary["message"])
        self.assertEqual(pair.snp_status, SNP_STATUS_SKIPPED)

    def test_no_genomic_coordinates_is_skipped(self):
        var_info = SimpleNamespace(genomic_pos=None)
        pair = _genomic_pair()
        summary = annotate_primer_pairs_with_snp_awareness(
            var_info, [pair], "GRCh37", enabled=True
        )
        self.assertEqual(summary["status"], SNP_STATUS_SKIPPED)
        self.assertIn("No genomic coordinates available", summary["message"])

    @patch("primer_designer_app.utils.snp_awareness.EnsemblClient")
    def test_ensembl_overlap_error_marks_pairs_error(self, mock_client_cls):
        mock_client_cls.return_value.get_overlapping_variations_for_region.side_effect = RuntimeError(
            "Ensembl unreachable"
        )
        pair = _genomic_pair()
        summary = annotate_primer_pairs_with_snp_awareness(
            _genomic_var_info(), [pair], "GRCh37", enabled=True
        )
        self.assertEqual(summary["status"], SNP_STATUS_ERROR)
        self.assertIn("Ensembl SNP lookup failed", summary["message"])
        self.assertEqual(pair.snp_status, SNP_STATUS_ERROR)
        self.assertEqual(pair.snp_conflicts, [])

    @patch("primer_designer_app.utils.snp_awareness.EnsemblClient")
    def test_variation_details_error_marks_pairs_error(self, mock_client_cls):
        mock_client_cls.return_value.get_overlapping_variations_for_region.return_value = [
            {"id": "rs100", "start": 1905, "end": 1905, "alleles": ["A", "G"]}
        ]
        mock_client_cls.return_value.get_variation_details_batch.side_effect = (
            RuntimeError("Ensembl frequency lookup unreachable")
        )
        pair = _genomic_pair()
        summary = annotate_primer_pairs_with_snp_awareness(
            _genomic_var_info(), [pair], "GRCh37", enabled=True
        )
        self.assertEqual(summary["status"], SNP_STATUS_ERROR)
        self.assertIn("Ensembl SNP frequency lookup failed", summary["message"])
        self.assertEqual(pair.snp_status, SNP_STATUS_ERROR)

    @patch("primer_designer_app.utils.snp_awareness.EnsemblClient")
    def test_success_with_binding_conflict(self, mock_client_cls):
        # Hit at genomic 1905 -> template index 1905 - region_start(1000) = 905,
        # which overlaps the forward primer window [900, 919].
        mock_client_cls.return_value.get_overlapping_variations_for_region.return_value = [
            {
                "id": "rs100",
                "start": 1905,
                "end": 1905,
                "alleles": ["A", "G"],
                "source": "dbSNP",
            }
        ]
        mock_client_cls.return_value.get_variation_details_batch.return_value = {
            "rs100": {"MAF": 0.05}
        }
        pair = _genomic_pair()
        summary = annotate_primer_pairs_with_snp_awareness(
            _genomic_var_info(), [pair], "GRCh37", enabled=True
        )
        self.assertEqual(summary["status"], "ok")
        self.assertEqual(summary["variant_count"], 1)
        self.assertEqual(summary["pairs_with_binding_conflict"], 1)
        self.assertEqual(pair.snp_status, SNP_STATUS_CONFLICT)
        self.assertEqual(len(pair.snp_conflicts), 1)
        self.assertEqual(pair.snp_conflicts[0]["primer"], "forward")

    @patch("primer_designer_app.utils.snp_awareness.EnsemblClient")
    def test_success_with_no_hits_after_maf_filter(self, mock_client_cls):
        mock_client_cls.return_value.get_overlapping_variations_for_region.return_value = [
            {"id": "rs_rare", "start": 1905, "end": 1905, "alleles": ["A", "G"]}
        ]
        mock_client_cls.return_value.get_variation_details_batch.return_value = {
            "rs_rare": {"MAF": 0.001}
        }
        pair = _genomic_pair()
        summary = annotate_primer_pairs_with_snp_awareness(
            _genomic_var_info(), [pair], "GRCh37", enabled=True
        )
        self.assertEqual(summary["status"], "ok")
        self.assertEqual(summary["variant_count"], 0)
        self.assertEqual(pair.snp_status, SNP_STATUS_NONE)


if __name__ == "__main__":
    unittest.main()
