import unittest

from primer_designer_app.exceptions import (
    InvalidReferenceSequenceError,
    NoPrimerPairsFoundError,
)
from primer_designer_app.utils.design_validation import (
    masked_fraction,
    validate_primer_search_results,
    validate_reference_sequence_for_design,
)
from primer_designer_app.utils.primer_utils import PrimerSearchResults


class DesignValidationTests(unittest.TestCase):
    def test_masked_fraction_all_n(self):
        self.assertEqual(masked_fraction("N" * 100), 1.0)

    def test_masked_fraction_mixed(self):
        seq = "A" * 50 + "N" * 50
        self.assertEqual(masked_fraction(seq), 0.5)

    def test_reject_mostly_masked_sequence(self):
        with self.assertRaises(InvalidReferenceSequenceError):
            validate_reference_sequence_for_design("N" * 2001)

    def test_accept_normal_sequence(self):
        validate_reference_sequence_for_design("ACGT" * 250)

    def test_reject_empty_primer_results(self):
        with self.assertRaises(NoPrimerPairsFoundError):
            validate_primer_search_results(PrimerSearchResults())

    def test_accept_nonempty_primer_results(self):
        from primer_designer_app.utils.primer_utils import PrimerPairResult

        results = PrimerSearchResults()
        results.primer_pairs = [
            PrimerPairResult(
                index=0,
                left_seq="A" * 20,
                right_seq="T" * 20,
                penalty=0.1,
                product_size=100,
            )
        ]
        validate_primer_search_results(results)
