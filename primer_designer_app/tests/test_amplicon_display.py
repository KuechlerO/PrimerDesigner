import unittest

from primer_designer_app.utils.amplicon_display import (
    amplicon_chrom_label,
    extract_amplicon_summary,
    format_penalty_score,
    truncate_product_seq,
)


class ExtractAmpliconSummaryTests(unittest.TestCase):
    def test_genomic_amplicon_uses_period_formatted_positions(self):
        amplicon = {"Chrom": "7", "ForPos": 19997582, "RevEnd": 19998000}
        self.assertEqual(
            extract_amplicon_summary(amplicon), "7: 19.997.582 - 19.998.000"
        )

    def test_transcriptomic_amplicon_uses_gene_and_transcript(self):
        # parts[-4] is the "gene symbol" slot per the pipe-delimited Chrom format.
        chrom = "ENST00000456328.2|GENE1|+|BRCA1|protein_coding"
        parts = chrom.split("|")
        amplicon = {"Chrom": chrom, "ForPos": 100, "RevEnd": 200}
        expected_gene_symbol = parts[-4]
        expected = f"{expected_gene_symbol} (ENST00000456328.2): 100 - 200"
        self.assertEqual(extract_amplicon_summary(amplicon), expected)

    def test_empty_chrom_uses_genomic_style_summary(self):
        amplicon = {"Chrom": "", "ForPos": 10, "RevEnd": 20}
        self.assertEqual(extract_amplicon_summary(amplicon), ": 10 - 20")


class AmpliconChromLabelTests(unittest.TestCase):
    def test_genomic_label_is_just_the_chromosome(self):
        self.assertEqual(amplicon_chrom_label({"Chrom": "X"}), "X")

    def test_transcriptomic_label_uses_gene_and_transcript(self):
        chrom = "ENST00000456328.2|GENE1|+|BRCA1|protein_coding"
        parts = chrom.split("|")
        expected_gene_symbol = parts[-4]
        self.assertEqual(
            amplicon_chrom_label({"Chrom": chrom}),
            f"{expected_gene_symbol} (ENST00000456328.2)",
        )

    def test_missing_chrom_returns_empty_string(self):
        self.assertEqual(amplicon_chrom_label({}), "")


class FormatPenaltyScoreTests(unittest.TestCase):
    def test_formats_to_two_decimals(self):
        self.assertEqual(format_penalty_score(1.5), "1.50")
        self.assertEqual(format_penalty_score(0), "0.00")

    def test_none_or_empty_returns_empty_string(self):
        self.assertEqual(format_penalty_score(None), "")
        self.assertEqual(format_penalty_score(""), "")

    def test_non_numeric_falls_back_to_str(self):
        self.assertEqual(format_penalty_score("N/A"), "N/A")


class TruncateProductSeqTests(unittest.TestCase):
    def test_short_sequence_unchanged(self):
        self.assertEqual(truncate_product_seq("ACGT", max_len=64), "ACGT")

    def test_long_sequence_truncated_with_ellipsis(self):
        seq = "A" * 100
        result = truncate_product_seq(seq, max_len=64)
        self.assertEqual(len(result), 64)
        self.assertTrue(result.endswith("\u2026"))
        self.assertEqual(result[:-1], "A" * 63)

    def test_none_returns_empty_string(self):
        self.assertEqual(truncate_product_seq(None), "")

    def test_default_max_len(self):
        seq = "A" * 200
        result = truncate_product_seq(seq)
        self.assertEqual(len(result), 64)


if __name__ == "__main__":
    unittest.main()
