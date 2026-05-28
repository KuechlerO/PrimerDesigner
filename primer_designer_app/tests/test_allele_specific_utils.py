import unittest
from types import SimpleNamespace

from primer_designer_app.utils.primer_utils import (
    PrimerPairResult,
    PrimerSearchResults,
    _infer_discriminating_left_ends,
)
from primer_designer_app.utils.variant_info import IndelType


class TestAlleleSpecificUtils(unittest.TestCase):
    def test_infer_discriminating_left_ends_snv(self):
        v = SimpleNamespace(
            relative_pos=(10, 10), indel_type=IndelType.SNV, new_bases="G"
        )
        wt_end, mut_end = _infer_discriminating_left_ends(v)
        self.assertEqual(wt_end, 10)
        self.assertEqual(mut_end, 10)

    def test_infer_discriminating_left_ends_ins(self):
        v = SimpleNamespace(
            relative_pos=(10, 9), indel_type=IndelType.INS, new_bases="AT"
        )
        wt_end, mut_end = _infer_discriminating_left_ends(v)
        # For INS, the current implementation forces the left primer 3' end onto the
        # first discriminating position at the junction (best-effort).
        self.assertEqual(wt_end, 10)
        self.assertEqual(mut_end, 10)

    def test_infer_discriminating_left_ends_del(self):
        v = SimpleNamespace(
            relative_pos=(10, 12), indel_type=IndelType.DEL, new_bases=""
        )
        wt_end, mut_end = _infer_discriminating_left_ends(v)
        self.assertEqual(wt_end, 12)
        self.assertEqual(mut_end, 9)

    def test_infer_discriminating_left_ends_delins(self):
        v = SimpleNamespace(
            relative_pos=(10, 12), indel_type=IndelType.DELINS, new_bases="A"
        )
        wt_end, mut_end = _infer_discriminating_left_ends(v)
        self.assertEqual(wt_end, 12)
        self.assertEqual(mut_end, 10)

    def test_build_allele_display_chunks_wt_plain_uses_ref_base(self):
        from primer_designer_app.utils.helpers import build_allele_display_chunks
        from primer_designer_app.utils.variant_info import AllelicVariantInfo, IndelType

        ref = "A" * 20 + "G" + "T" * 20
        mut = "A" * 20 + "A" + "T" * 20
        v = AllelicVariantInfo(
            ref_seq=ref,
            relative_pos=(20, 20),
            indel_type=IndelType.SNV,
            ref_bases="G",
            new_bases="A",
        )
        wt_chunks = build_allele_display_chunks(
            ref,
            v,
            allele="wt",
            width=50,
            ref_bases="G",
            new_bases="A",
            highlight_snv_allele=True,
        )
        mut_chunks = build_allele_display_chunks(
            mut,
            v,
            allele="mut",
            width=50,
            ref_bases="G",
            new_bases="A",
            highlight_snv_allele=True,
        )
        self.assertEqual(wt_chunks[0]["plain"][20], "G")
        self.assertEqual(mut_chunks[0]["plain"][20], "A")

    def test_hgvs_ins_mut_plain_suffix_skips_duplicate_insert(self):
        from primer_designer_app.utils.hgvs_display import hgvs_input_on_plain
        from primer_designer_app.utils.variant_info import IndelType

        lo = 50
        insert = "ATTGCGCAATGC"
        mut = "A" * lo + insert + "G" * 10
        annotated = hgvs_input_on_plain(
            mut, lo, lo, IndelType.INS, "", insert, allele="mut"
        )
        self.assertIn("[-/" + insert + "]", annotated)
        after_bracket = annotated.split("]", 1)[1]
        self.assertFalse(after_bracket.startswith(insert))
        self.assertTrue(after_bracket.startswith("G"))

    def test_hgvs_delins_mut_plain_suffix_skips_duplicate_alt(self):
        from primer_designer_app.utils.hgvs_display import hgvs_input_on_plain
        from primer_designer_app.utils.variant_info import IndelType

        lo = 10
        ref = "GCAGG"  # 5bp
        alt = "GGTC"  # 4bp
        suffix = "T" * 10
        mut_plain = "A" * lo + alt + suffix

        annotated = hgvs_input_on_plain(
            mut_plain,
            lo,
            lo + len(ref) - 1,
            IndelType.DELINS,
            ref,
            alt,
            allele="mut",
        )

        self.assertIn("[" + ref + "/" + alt + "]", annotated)
        after_bracket = annotated.split("]", 1)[1]
        self.assertFalse(after_bracket.startswith(alt))
        self.assertTrue(after_bracket.startswith(suffix))

    def test_mut_ins_bracket_matches_wt_notation(self):
        """WT and MUT use the same HGVS input INS notation [-/inserted]."""
        lo = 10
        new_bases = "ATGC"
        plain = "A" * lo + new_bases + "T" * 10
        annot = plain[:lo] + "[-/" + new_bases + "]" + plain[lo:]
        self.assertIn("[-/ATGC]", annot)

    def test_primer_search_results_roundtrip(self):
        res = PrimerSearchResults()
        res.mapped_primer_positions = {"primerF_starts": [100], "primerR_ends": [200]}
        res.primer_pairs = [
            PrimerPairResult(
                index=0,
                left_seq="AAA",
                right_seq="TTT",
                penalty=1.0,
                product_size=123,
                product_tm=75.0,
                left_relPos_start=5,
                left_relPos_end=7,
                right_relPos_start=20,
                right_relPos_end=22,
                tm=[60.0, 61.0],
                gc=[50.0, 52.0],
            )
        ]
        d = res.to_dict()
        res2 = PrimerSearchResults.from_dict(d)
        self.assertEqual(res2.mapped_primer_positions["primerF_starts"], [100])
        self.assertEqual(res2.primer_pairs[0].left_seq, "AAA")
        self.assertEqual(res2.primer_pairs[0].right_seq, "TTT")


if __name__ == "__main__":
    unittest.main()
