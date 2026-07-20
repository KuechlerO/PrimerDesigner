import unittest
from dataclasses import dataclass

from primer_designer_app.utils.primer3_post import (
    PRIMER3_OVERRIDE_FIELDS,
    _coerce,
    parse_primer3_overrides_from_post,
)


@dataclass
class _FakeRequest:
    POST: dict


class CoerceTests(unittest.TestCase):
    def test_coerce_int(self):
        self.assertEqual(_coerce("18", "int"), 18)
        self.assertIsInstance(_coerce("18", "int"), int)

    def test_coerce_float(self):
        self.assertEqual(_coerce("60.5", "float"), 60.5)
        self.assertIsInstance(_coerce("60.5", "float"), float)

    def test_coerce_strips_whitespace(self):
        self.assertEqual(_coerce("  20  ", "int"), 20)

    def test_coerce_invalid_int_raises(self):
        with self.assertRaises(ValueError):
            _coerce("abc", "int")

    def test_coerce_invalid_float_raises(self):
        with self.assertRaises(ValueError):
            _coerce("abc", "float")

    def test_coerce_unknown_kind_raises(self):
        with self.assertRaises(ValueError):
            _coerce("20", "bogus_kind")


class ParsePrimer3OverridesFromPostTests(unittest.TestCase):
    def test_empty_fields_are_skipped(self):
        request = _FakeRequest(
            POST={"p3_PRIMER_MIN_SIZE": "   ", "p3_PRIMER_MAX_SIZE": ""}
        )
        overrides = parse_primer3_overrides_from_post(request)
        self.assertEqual(overrides, {})

    def test_missing_fields_are_skipped(self):
        request = _FakeRequest(POST={})
        overrides = parse_primer3_overrides_from_post(request)
        self.assertEqual(overrides, {})

    def test_invalid_values_are_skipped(self):
        request = _FakeRequest(
            POST={"p3_PRIMER_MIN_SIZE": "not-a-number", "p3_PRIMER_MIN_TM": "58.5"}
        )
        overrides = parse_primer3_overrides_from_post(request)
        self.assertNotIn("PRIMER_MIN_SIZE", overrides)
        self.assertEqual(overrides["PRIMER_MIN_TM"], 58.5)

    def test_valid_values_are_coerced_to_correct_types(self):
        request = _FakeRequest(
            POST={
                "p3_PRIMER_MIN_SIZE": "18",
                "p3_PRIMER_MAX_SIZE": "22",
                "p3_PRIMER_MIN_TM": "58.0",
                "p3_PRIMER_MAX_TM": "62.0",
                "p3_PRIMER_GC_CLAMP": "1",
            }
        )
        overrides = parse_primer3_overrides_from_post(request)
        self.assertEqual(overrides["PRIMER_MIN_SIZE"], 18)
        self.assertIsInstance(overrides["PRIMER_MIN_SIZE"], int)
        self.assertEqual(overrides["PRIMER_MAX_SIZE"], 22)
        self.assertEqual(overrides["PRIMER_MIN_TM"], 58.0)
        self.assertIsInstance(overrides["PRIMER_MIN_TM"], float)
        self.assertEqual(overrides["PRIMER_GC_CLAMP"], 1)

    def test_only_whitelisted_keys_are_considered(self):
        request = _FakeRequest(POST={"p3_NOT_A_REAL_KEY": "123"})
        overrides = parse_primer3_overrides_from_post(request)
        self.assertEqual(overrides, {})

    def test_all_whitelisted_fields_round_trip(self):
        post = {f"p3_{key}": "1" for key, _kind in PRIMER3_OVERRIDE_FIELDS}
        request = _FakeRequest(POST=post)
        overrides = parse_primer3_overrides_from_post(request)
        self.assertEqual(len(overrides), len(PRIMER3_OVERRIDE_FIELDS))
        for key, kind in PRIMER3_OVERRIDE_FIELDS:
            expected_type = int if kind == "int" else float
            self.assertIsInstance(overrides[key], expected_type)


if __name__ == "__main__":
    unittest.main()
