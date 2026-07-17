import unittest
from types import SimpleNamespace

from primer_designer_app.utils.primer_utils import (
    INSILICO_NOT_APPLICABLE,
    INSILICO_OK,
    INSILICO_OK_EMPTY,
    PrimerPairResult,
    _infer_legacy_insilico_status,
    _pick_common_reverse_from_primer3,
    build_primer3_global_args,
    get_primers_from_primer3,
    primer_pair_from_dict,
)


class GetPrimersFromPrimer3Tests(unittest.TestCase):
    def _raw_pair(self, i: int, *, left_start, left_len, right_5p, right_len):
        return {
            f"PRIMER_LEFT_{i}": [left_start, left_len],
            f"PRIMER_LEFT_{i}_SEQUENCE": "A" * left_len,
            f"PRIMER_RIGHT_{i}": [right_5p, right_len],
            f"PRIMER_RIGHT_{i}_SEQUENCE": "T" * right_len,
            f"PRIMER_LEFT_{i}_GC_PERCENT": 50.0,
            f"PRIMER_RIGHT_{i}_GC_PERCENT": 55.0,
            f"PRIMER_LEFT_{i}_TM": 60.0,
            f"PRIMER_RIGHT_{i}_TM": 61.0,
            f"PRIMER_PAIR_{i}_PENALTY": 0.5 + i,
            f"PRIMER_PAIR_{i}_PRODUCT_SIZE": 100 + i,
            f"PRIMER_PAIR_{i}_PRODUCT_TM": 75.0,
        }

    def test_well_formed_two_pairs_right_primer_0based_coords(self):
        dicey_primer = {"PRIMER_PAIR_NUM_RETURNED": 2}
        dicey_primer.update(
            self._raw_pair(0, left_start=10, left_len=20, right_5p=100, right_len=20)
        )
        dicey_primer.update(
            self._raw_pair(1, left_start=15, left_len=18, right_5p=200, right_len=22)
        )

        pairs = get_primers_from_primer3(dicey_primer)

        self.assertEqual(len(pairs), 2)

        pair0 = pairs[0]
        self.assertEqual(pair0.left_relPos_start, 10)
        self.assertEqual(pair0.left_relPos_end, 29)
        # right primer: 5' position 100 (0-based), length 20 -> start=81, end=100
        self.assertEqual(pair0.right_relPos_start, 81)
        self.assertEqual(pair0.right_relPos_end, 100)

        pair1 = pairs[1]
        self.assertEqual(pair1.left_relPos_start, 15)
        self.assertEqual(pair1.left_relPos_end, 32)
        # 5' position 200, length 22 -> start = 200 - 21 = 179, end = 200
        self.assertEqual(pair1.right_relPos_start, 179)
        self.assertEqual(pair1.right_relPos_end, 200)

    def test_zero_pairs_returned(self):
        dicey_primer = {"PRIMER_PAIR_NUM_RETURNED": 0}
        pairs = get_primers_from_primer3(dicey_primer)
        self.assertEqual(pairs, [])

    def test_missing_num_returned_defaults_to_zero(self):
        self.assertEqual(get_primers_from_primer3({}), [])

    def test_skips_malformed_pair(self):
        dicey_primer = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0": ["not-an-int", 20],
            "PRIMER_LEFT_0_SEQUENCE": "AAAA",
            "PRIMER_RIGHT_0": [100, 20],
            "PRIMER_RIGHT_0_SEQUENCE": "TTTT",
        }
        pairs = get_primers_from_primer3(dicey_primer)
        self.assertEqual(pairs, [])


class BuildPrimer3GlobalArgsTests(unittest.TestCase):
    def _prim_set(self, **overrides):
        return SimpleNamespace(
            tm=60,
            gc=50,
            max_poly_x=4,
            productsize_range=[400, 800],
            primer3_overrides=overrides,
        )

    def test_base_defaults_without_overrides(self):
        prim_set = self._prim_set()
        args = build_primer3_global_args(prim_set)
        self.assertEqual(args["PRIMER_OPT_TM"], 60.0)
        self.assertEqual(args["PRIMER_MIN_TM"], 58.0)
        self.assertEqual(args["PRIMER_MAX_TM"], 62.0)
        self.assertEqual(args["PRIMER_MIN_SIZE"], 18)
        self.assertEqual(args["PRIMER_PRODUCT_SIZE_RANGE"], [400, 800])

    def test_whitelisted_override_merges_and_wins(self):
        prim_set = self._prim_set(PRIMER_MIN_SIZE=25)
        args = build_primer3_global_args(prim_set)
        self.assertEqual(args["PRIMER_MIN_SIZE"], 25)
        # Untouched base defaults remain.
        self.assertEqual(args["PRIMER_MAX_SIZE"], 22)

    def test_non_whitelisted_override_is_rejected(self):
        prim_set = self._prim_set(PRIMER_MIN_SIZE=25, NOT_A_REAL_KEY=999)
        args = build_primer3_global_args(prim_set)
        self.assertEqual(args["PRIMER_MIN_SIZE"], 25)
        self.assertNotIn("NOT_A_REAL_KEY", args)

    def test_missing_primer3_overrides_attribute_is_ok(self):
        prim_set = SimpleNamespace(
            tm=60, gc=50, max_poly_x=4, productsize_range=[400, 800]
        )
        args = build_primer3_global_args(prim_set)
        self.assertEqual(args["PRIMER_OPT_TM"], 60.0)


class InferLegacyInsilicoStatusTests(unittest.TestCase):
    def test_none_amplicons_is_ok_empty(self):
        self.assertEqual(_infer_legacy_insilico_status({}), INSILICO_OK_EMPTY)

    def test_na_marker_is_not_applicable(self):
        self.assertEqual(
            _infer_legacy_insilico_status({"amplicons": ["N/A"]}),
            INSILICO_NOT_APPLICABLE,
        )

    def test_dict_amplicons_is_ok(self):
        self.assertEqual(
            _infer_legacy_insilico_status({"amplicons": [{"Seq": "ACGT"}]}),
            INSILICO_OK,
        )

    def test_empty_list_is_ok_empty(self):
        self.assertEqual(
            _infer_legacy_insilico_status({"amplicons": []}), INSILICO_OK_EMPTY
        )


class PrimerPairFromDictTests(unittest.TestCase):
    def test_infers_missing_insilico_status(self):
        pair = primer_pair_from_dict(
            {
                "index": 0,
                "left_seq": "AAA",
                "right_seq": "TTT",
                "penalty": 0.1,
                "product_size": 100,
                "amplicons": [{"Seq": "ACGT"}],
            }
        )
        self.assertIsInstance(pair, PrimerPairResult)
        self.assertEqual(pair.insilico_status, INSILICO_OK)
        self.assertIsNone(pair.insilico_error_detail)
        self.assertIsNone(pair.snp_status)
        self.assertIsNone(pair.snp_conflicts)

    def test_preserves_explicit_insilico_status(self):
        pair = primer_pair_from_dict(
            {
                "index": 0,
                "left_seq": "AAA",
                "right_seq": "TTT",
                "penalty": 0.1,
                "product_size": 100,
                "insilico_status": "error",
            }
        )
        self.assertEqual(pair.insilico_status, "error")


class PickCommonReverseFromPrimer3Tests(unittest.TestCase):
    def test_no_candidates_returns_none(self):
        raw = {"PRIMER_RIGHT_NUM_RETURNED": 0}
        self.assertIsNone(
            _pick_common_reverse_from_primer3(raw, min_left_end=100, min_gap=80)
        )

    def test_candidate_below_min_start_is_skipped(self):
        raw = {
            "PRIMER_RIGHT_NUM_RETURNED": 1,
            "PRIMER_RIGHT_0": [150, 20],
            "PRIMER_RIGHT_0_SEQUENCE": "T" * 20,
        }
        # min_start = 100 + 80 = 180; candidate start (150) is not > 180.
        self.assertIsNone(
            _pick_common_reverse_from_primer3(raw, min_left_end=100, min_gap=80)
        )

    def test_valid_candidate_is_returned(self):
        raw = {
            "PRIMER_RIGHT_NUM_RETURNED": 2,
            "PRIMER_RIGHT_0": [150, 20],
            "PRIMER_RIGHT_0_SEQUENCE": "T" * 20,
            "PRIMER_RIGHT_1": [300, 20],
            "PRIMER_RIGHT_1_SEQUENCE": "G" * 20,
        }
        result = _pick_common_reverse_from_primer3(raw, min_left_end=100, min_gap=80)
        self.assertEqual(result, ("G" * 20, 300))

    def test_missing_position_or_sequence_is_skipped(self):
        raw = {
            "PRIMER_RIGHT_NUM_RETURNED": 1,
            "PRIMER_RIGHT_0": None,
            "PRIMER_RIGHT_0_SEQUENCE": None,
        }
        self.assertIsNone(
            _pick_common_reverse_from_primer3(raw, min_left_end=100, min_gap=80)
        )


if __name__ == "__main__":
    unittest.main()
