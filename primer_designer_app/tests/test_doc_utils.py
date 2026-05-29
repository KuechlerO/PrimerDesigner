import unittest
from dataclasses import replace
from types import SimpleNamespace

from primer_designer_app.utils.display_utils import (
    REPORT_DISPLAY_FLANK,
    compute_report_display_bounds,
)
from primer_designer_app.utils.doc_utils import (
    HIGHLIGHT_PRIMER,
    HIGHLIGHT_SNP,
    HIGHLIGHT_SNP_CONFLICT,
    HIGHLIGHT_VCF,
    build_template_highlight_lookup,
    prepare_report_sequence_view,
)
from primer_designer_app.utils.variant_info import IndelType
from primer_designer_app.utils.primer_utils import PrimerPairResult


class TestDocUtilsHighlights(unittest.TestCase):
    def _pair(self) -> PrimerPairResult:
        return PrimerPairResult(
            index=0,
            left_seq="AAA",
            right_seq="TTT",
            penalty=1.0,
            product_size=100,
            product_tm=70.0,
            left_relPos_start=10,
            left_relPos_end=19,
            right_relPos_start=80,
            right_relPos_end=89,
            tm=[60.0, 61.0],
            gc=[50.0, 52.0],
        )

    def test_vcf_and_snp_layers(self):
        lookup = build_template_highlight_lookup(
            100,
            vcf_hits=[{"template_start": 5, "template_end": 5}],
            snp_hits=[{"template_start": 50, "template_end": 50}],
            primer_pair=self._pair(),
        )
        self.assertEqual(lookup[5], HIGHLIGHT_VCF)
        self.assertEqual(lookup[50], HIGHLIGHT_SNP)
        self.assertEqual(lookup[15], HIGHLIGHT_PRIMER)

    def test_snp_conflict_uses_orange(self):
        lookup = build_template_highlight_lookup(
            100,
            vcf_hits=[],
            snp_hits=[{"template_start": 12, "template_end": 12}],
            primer_pair=self._pair(),
        )
        self.assertEqual(lookup[12], "orange")
        self.assertEqual(lookup[12], HIGHLIGHT_SNP_CONFLICT)

    def test_skip_user_variant_interval_for_snp(self):
        lookup = build_template_highlight_lookup(
            100,
            vcf_hits=[],
            snp_hits=[{"template_start": 20, "template_end": 20}],
            primer_pair=self._pair(),
            skip_snp_interval=(20, 24),
        )
        self.assertIsNone(lookup[20])

    def test_prepare_report_sequence_view_slices_to_region(self):
        plain = "N" * 2000
        var = SimpleNamespace(
            relative_pos=(1000, 1000),
            indel_type=IndelType.SNV,
            ref_bases="A",
            new_bases="G",
        )
        settings = SimpleNamespace(target=(995, 10))
        pair = replace(
            self._pair(),
            left_relPos_start=980,
            left_relPos_end=989,
            right_relPos_start=1010,
            right_relPos_end=1019,
        )
        _, plain_slice, shifted, _, _, _ = prepare_report_sequence_view(
            var,
            settings,
            pair,
            design_template=plain,
            allele="wt",
            vcf_hits=[],
            snp_hits=[],
        )
        self.assertLess(len(plain_slice), len(plain))
        self.assertGreater(len(plain_slice), 0)
        self.assertEqual(shifted.left_relPos_start, 980 - (980 - REPORT_DISPLAY_FLANK))
        self.assertGreaterEqual(shifted.right_relPos_end, 0)

    def test_prepare_report_delins_uses_mut_allele_suffix(self):
        from primer_designer_app.utils.hgvs_display import hgvs_input_on_plain

        lo, ref, alt = 100, "GCAGG", "GGTC"
        plain = "A" * lo + alt + "ACGTACGTACGT"
        hi = lo + len(ref) - 1
        annot_wt = hgvs_input_on_plain(
            plain, lo, hi, IndelType.DELINS, ref, alt, allele="wt"
        )
        annot_mut = hgvs_input_on_plain(
            plain, lo, hi, IndelType.DELINS, ref, alt, allele="mut"
        )
        suffix_plain = plain[lo + len(alt) : lo + len(alt) + 6]
        self.assertEqual(annot_mut.split("]", 1)[1][:6], suffix_plain)
        self.assertNotEqual(annot_wt.split("]", 1)[1][:6], suffix_plain)

    def test_report_window_spans_both_primers(self):
        pair = replace(
            self._pair(),
            left_relPos_start=80,
            left_relPos_end=99,
            right_relPos_start=700,
            right_relPos_end=719,
            product_size=650,
        )
        start, end = compute_report_display_bounds(2000, 1000, 1000, 950, 101, pair)
        self.assertLessEqual(start, 80)
        self.assertGreaterEqual(end, 720)
        self.assertGreater(end - start, 700)


if __name__ == "__main__":
    unittest.main()
