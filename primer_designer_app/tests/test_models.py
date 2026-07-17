from django.test import TestCase

from primer_designer_app.models import DesignResultsSummary, PrimerSettingsModel
from primer_designer_app.utils.primer_utils import PrimerPairResult, PrimerSearchResults
from primer_designer_app.utils.variant_info import (
    IndelType,
    SequenceVariantInfo,
    StructuralVariantInfo,
)


def _make_primer_settings(**overrides) -> PrimerSettingsModel:
    defaults = dict(
        tm=60,
        gc=50,
        max_poly_x=4,
        productsize_range=[400, 800],
        reference_genome="GRCh37",
    )
    defaults.update(overrides)
    return PrimerSettingsModel.objects.create(**defaults)


class PrimerSettingsModelSetTargetTests(TestCase):
    def test_set_target_computes_padded_window_and_saves(self):
        primer_settings = _make_primer_settings(target_padding=50)
        primer_settings.set_target((100, 150))

        self.assertEqual(primer_settings.target, [50, 150])

        reloaded = PrimerSettingsModel.objects.get(pk=primer_settings.pk)
        self.assertEqual(reloaded.target, [50, 150])

    def test_set_target_rejects_padding_below_one(self):
        primer_settings = _make_primer_settings(target_padding=0)
        with self.assertRaisesRegex(ValueError, "Invalid target_padding"):
            primer_settings.set_target((100, 150))

    def test_set_target_rejects_padding_above_500(self):
        primer_settings = _make_primer_settings(target_padding=501)
        with self.assertRaisesRegex(ValueError, "Invalid target_padding"):
            primer_settings.set_target((100, 150))


class PrimerSettingsModelSetContextTests(TestCase):
    def test_set_context_valid_values_are_saved(self):
        primer_settings = _make_primer_settings()
        primer_settings.set_context("transcriptomic")

        self.assertEqual(primer_settings.context, "transcriptomic")
        reloaded = PrimerSettingsModel.objects.get(pk=primer_settings.pk)
        self.assertEqual(reloaded.context, "transcriptomic")

    def test_set_context_invalid_value_raises(self):
        primer_settings = _make_primer_settings()
        with self.assertRaisesRegex(ValueError, "Invalid context"):
            primer_settings.set_context("bogus")
        # Unchanged on failure.
        self.assertEqual(primer_settings.context, "genomic")


class DesignResultsSummarySequenceVariantRoundtripTests(TestCase):
    def test_save_and_get_variant_info_for_sequence_variant(self):
        var_info = SequenceVariantInfo(input_seq="AC[G>A]GT", ref_genome="GRCh37")
        primer_settings = _make_primer_settings()

        results = PrimerSearchResults()
        results.primer_pairs = [
            PrimerPairResult(
                index=0,
                left_seq="A" * 20,
                right_seq="T" * 20,
                penalty=0.1,
                product_size=100,
                left_relPos_start=0,
                left_relPos_end=19,
                right_relPos_start=80,
                right_relPos_end=99,
            )
        ]

        summary = DesignResultsSummary()
        summary.save_primer_results(
            results, primer_settings, var_info, snp_analysis_data={"status": "ok"}
        )

        self.assertIsNotNone(summary.pk)
        reloaded = DesignResultsSummary.objects.get(pk=summary.pk)

        self.assertFalse(reloaded.is_structural_variant_design())

        restored_var_info = reloaded.get_variant_info()
        self.assertIsInstance(restored_var_info, SequenceVariantInfo)
        self.assertEqual(restored_var_info.ref_seq, "ACGGT")
        self.assertEqual(restored_var_info.indel_type, IndelType.SNV)
        self.assertEqual(restored_var_info.relative_pos, [2, 2])

        restored_results = reloaded.get_primer_search_results()
        self.assertIsInstance(restored_results, PrimerSearchResults)
        self.assertEqual(len(restored_results.primer_pairs), 1)
        self.assertEqual(restored_results.primer_pairs[0].left_seq, "A" * 20)

        self.assertEqual(reloaded.snp_analysis_data, {"status": "ok"})


class DesignResultsSummaryStructuralVariantRoundtripTests(TestCase):
    def test_save_and_get_structural_variant_results(self):
        sv_info = StructuralVariantInfo(
            chromosome="1",
            start_position=1000,
            end_position=5000,
            reference_genome="GRCh37",
        )
        windows = sv_info.create_design_windows()

        sv_results = {}
        for i, window in enumerate(windows):
            pair = PrimerPairResult(
                index=0,
                left_seq="A" * 20,
                right_seq="T" * 20,
                penalty=0.1 + i,
                product_size=100,
                left_relPos_start=0,
                left_relPos_end=19,
                right_relPos_start=80,
                right_relPos_end=99,
            )
            sv_results[window.label] = {
                "design_window": window,
                "primer_rows": [
                    {
                        "pair": pair,
                        "genomic_positions": {
                            "forward_start": window.window_start_genomic,
                            "forward_end": window.window_start_genomic + 19,
                            "reverse_start": window.window_start_genomic + 80,
                            "reverse_end": window.window_start_genomic + 99,
                        },
                    }
                ],
            }

        # save_structural_variant_results saves primer_settings itself; no need to
        # create it beforehand.
        primer_settings = PrimerSettingsModel(
            tm=60,
            gc=50,
            max_poly_x=4,
            productsize_range=[400, 800],
            reference_genome="GRCh37",
        )

        summary = DesignResultsSummary()
        summary.save_structural_variant_results(primer_settings, sv_info, sv_results)

        self.assertIsNotNone(summary.pk)
        self.assertIsNotNone(primer_settings.pk)

        reloaded = DesignResultsSummary.objects.get(pk=summary.pk)
        self.assertTrue(reloaded.is_structural_variant_design())
        self.assertIsNone(reloaded.get_variant_info())
        self.assertIsNone(reloaded.get_primer_search_results())

        sv_data = reloaded.get_structural_variant_info_data()
        self.assertEqual(sv_data["chromosome"], "1")
        self.assertEqual(sv_data["start_position"], 1000)
        self.assertEqual(sv_data["end_position"], 5000)

        restored_results = reloaded.get_sv_primer_results()
        self.assertIn("upstream", restored_results)
        upstream_row = restored_results["upstream"]["primer_rows"][0]
        self.assertEqual(upstream_row["pair"].left_seq, "A" * 20)
        self.assertEqual(
            upstream_row["genomic_positions"]["forward_start"],
            sv_info.windows[0].window_start_genomic,
        )
