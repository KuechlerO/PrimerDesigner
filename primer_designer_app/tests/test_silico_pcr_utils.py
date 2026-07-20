import unittest
from unittest.mock import patch

from primer_designer_app.exceptions import InvalidInputError
from primer_designer_app.utils.primer_utils import INSILICO_OK
from primer_designer_app.utils.silico_pcr_utils import (
    MAX_PRIMER_LEN,
    MIN_PRIMER_LEN,
    build_primer_pair,
    build_silico_settings,
    normalize_primer_sequence,
    run_silico_pcr,
)


class NormalizePrimerSequenceTests(unittest.TestCase):
    def test_strips_whitespace_and_uppercases(self):
        self.assertEqual(
            normalize_primer_sequence(" atgc atgc at ", "Forward primer"),
            "ATGCATGCAT",
        )

    def test_rejects_empty(self):
        with self.assertRaises(InvalidInputError):
            normalize_primer_sequence("   ", "Forward primer")

    def test_rejects_invalid_alphabet(self):
        with self.assertRaises(InvalidInputError) as ctx:
            normalize_primer_sequence("ATGXATGCAT", "Reverse primer")
        self.assertIn("A, C, G, T, or N", str(ctx.exception))

    def test_rejects_too_short(self):
        with self.assertRaises(InvalidInputError) as ctx:
            normalize_primer_sequence("ATGC", "Forward primer")
        self.assertIn(str(MIN_PRIMER_LEN), str(ctx.exception))

    def test_rejects_too_long(self):
        seq = "A" * (MAX_PRIMER_LEN + 1)
        with self.assertRaises(InvalidInputError) as ctx:
            normalize_primer_sequence(seq, "Forward primer")
        self.assertIn(str(MAX_PRIMER_LEN), str(ctx.exception))

    def test_accepts_n_and_boundary_lengths(self):
        short = "A" * MIN_PRIMER_LEN
        long = "N" * MAX_PRIMER_LEN
        self.assertEqual(normalize_primer_sequence(short, "Forward primer"), short)
        self.assertEqual(normalize_primer_sequence(long, "Reverse primer"), long)


class BuildSilicoSettingsTests(unittest.TestCase):
    def test_genome_context(self):
        s = build_silico_settings("GRCh38", "genome")
        self.assertEqual(s.reference_genome, "GRCh38")
        self.assertEqual(s.context, "genomic")
        self.assertTrue(s.do_insilico_pcr)

    def test_transcriptome_context(self):
        s = build_silico_settings("GRCh37", "transcriptome")
        self.assertEqual(s.context, "transcriptomic")

    def test_rejects_none_context(self):
        with self.assertRaises(InvalidInputError):
            build_silico_settings("GRCh37", "none")

    def test_rejects_bad_genome(self):
        with self.assertRaises(InvalidInputError):
            build_silico_settings("hg19", "genome")


class RunSilicoPcrTests(unittest.TestCase):
    def test_builds_pair_and_runs_analysis(self):
        fwd = "ATGCGATCGATCGATCGATC"
        rev = "TACGTAGCTAGCTAGCTAGC"

        def fake_analysis(settings, pairs):
            pairs[0].amplicons = [
                {"Length": 120, "Chrom": "1", "ForPos": 10, "RevEnd": 130}
            ]
            pairs[0].insilico_status = INSILICO_OK
            pairs[0].insilico_error_detail = None

        with patch(
            "primer_designer_app.utils.silico_pcr_utils.do_insilico_analysis",
            side_effect=fake_analysis,
        ):
            pair, settings, note = run_silico_pcr(
                forward_raw=fwd,
                reverse_raw=rev,
                reference_genome="GRCh37",
                amplicon_check="genome",
            )

        self.assertEqual(pair.left_seq, fwd)
        self.assertEqual(pair.right_seq, rev)
        self.assertEqual(pair.insilico_status, INSILICO_OK)
        self.assertEqual(len(pair.amplicons), 1)
        self.assertEqual(settings.context, "genomic")
        self.assertIn("Genomic amplicons", note)

    def test_build_primer_pair_placeholders(self):
        pair = build_primer_pair("ATGCGATCGATCGATCGATC", "TACGTAGCTAGCTAGCTAGC")
        self.assertEqual(pair.penalty, 0.0)
        self.assertEqual(pair.product_size, 0)
        self.assertIsNone(pair.left_relPos_start)


if __name__ == "__main__":
    unittest.main()
