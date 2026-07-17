import unittest

from primer_designer_app.templatetags.custom_filters import (
    genomic_coord,
    genomic_range,
    insilico_cell_class,
    penalty_two_decimals,
    snp_cell_class,
    snp_conflict_summary,
    snp_status_label,
)
from primer_designer_app.utils.primer_utils import (
    INSILICO_ERROR,
    INSILICO_NOT_APPLICABLE,
    INSILICO_OK,
    INSILICO_OK_EMPTY,
)
from primer_designer_app.utils.snp_awareness import (
    SNP_STATUS_CAUTION,
    SNP_STATUS_CONFLICT,
    SNP_STATUS_ERROR,
    SNP_STATUS_NONE,
    SNP_STATUS_SKIPPED,
)


class GenomicCoordFilterTests(unittest.TestCase):
    def test_formats_with_periods(self):
        self.assertEqual(genomic_coord(19997582), "19.997.582")

    def test_none_returns_empty_string(self):
        self.assertEqual(genomic_coord(None), "")


class GenomicRangeFilterTests(unittest.TestCase):
    def test_formats_both_sides(self):
        self.assertEqual(genomic_range(1000000, 2000000), "1.000.000 - 2.000.000")

    def test_missing_end_returns_start_only(self):
        self.assertEqual(genomic_range(1000000, None), "1.000.000")


class SnpCellClassTests(unittest.TestCase):
    def test_known_statuses(self):
        self.assertEqual(snp_cell_class(SNP_STATUS_NONE), "snp-cell snp-cell--none")
        self.assertEqual(
            snp_cell_class(SNP_STATUS_CAUTION), "snp-cell snp-cell--caution"
        )
        self.assertEqual(
            snp_cell_class(SNP_STATUS_CONFLICT), "snp-cell snp-cell--conflict"
        )
        self.assertEqual(
            snp_cell_class(SNP_STATUS_SKIPPED), "snp-cell snp-cell--skipped"
        )
        self.assertEqual(snp_cell_class(SNP_STATUS_ERROR), "snp-cell snp-cell--error")

    def test_unknown_status_falls_back(self):
        self.assertEqual(snp_cell_class("bogus"), "snp-cell snp-cell--na")
        self.assertEqual(snp_cell_class(None), "snp-cell snp-cell--na")


class SnpStatusLabelTests(unittest.TestCase):
    def test_known_statuses(self):
        self.assertEqual(snp_status_label(SNP_STATUS_NONE), "Clear")
        self.assertEqual(snp_status_label(SNP_STATUS_CAUTION), "Caution")
        self.assertEqual(snp_status_label(SNP_STATUS_CONFLICT), "Conflict")
        self.assertEqual(snp_status_label(SNP_STATUS_SKIPPED), "N/A")
        self.assertEqual(snp_status_label(SNP_STATUS_ERROR), "Error")

    def test_unknown_status_falls_back(self):
        self.assertEqual(snp_status_label("bogus"), "\u2014")


class SnpConflictSummaryTests(unittest.TestCase):
    def test_empty_conflicts_returns_empty_string(self):
        self.assertEqual(snp_conflict_summary([]), "")
        self.assertEqual(snp_conflict_summary(None), "")

    def test_up_to_three_conflicts_lists_all(self):
        conflicts = [
            {"id": "rs1", "primer": "forward", "alleles": "A/G"},
            {"id": "rs2", "primer": "reverse", "alleles": "C/T"},
        ]
        summary = snp_conflict_summary(conflicts)
        self.assertIn("rs1 (forward, A/G)", summary)
        self.assertIn("rs2 (reverse, C/T)", summary)
        self.assertNotIn("more", summary)

    def test_more_than_three_conflicts_shows_overflow_count(self):
        conflicts = [
            {"id": f"rs{i}", "primer": "forward", "alleles": "A/G"} for i in range(5)
        ]
        summary = snp_conflict_summary(conflicts)
        # Only first 3 are listed explicitly.
        self.assertIn("rs0", summary)
        self.assertIn("rs1", summary)
        self.assertIn("rs2", summary)
        self.assertNotIn("rs3", summary)
        self.assertIn("(+2 more)", summary)


class InsilicoCellClassTests(unittest.TestCase):
    def test_known_statuses(self):
        self.assertEqual(
            insilico_cell_class(INSILICO_OK), "amplicon-cell amplicon-cell--ok"
        )
        self.assertEqual(
            insilico_cell_class(INSILICO_OK_EMPTY),
            "amplicon-cell amplicon-cell--empty",
        )
        self.assertEqual(
            insilico_cell_class(INSILICO_NOT_APPLICABLE),
            "amplicon-cell amplicon-cell--na",
        )
        self.assertEqual(
            insilico_cell_class(INSILICO_ERROR), "amplicon-cell amplicon-cell--error"
        )

    def test_unknown_status_falls_back(self):
        self.assertEqual(
            insilico_cell_class("bogus"), "amplicon-cell amplicon-cell--unknown"
        )
        self.assertEqual(
            insilico_cell_class(None), "amplicon-cell amplicon-cell--unknown"
        )


class PenaltyTwoDecimalsTests(unittest.TestCase):
    def test_formats_number(self):
        self.assertEqual(penalty_two_decimals(1.234), "1.23")

    def test_none_returns_empty_string(self):
        self.assertEqual(penalty_two_decimals(None), "")


if __name__ == "__main__":
    unittest.main()
