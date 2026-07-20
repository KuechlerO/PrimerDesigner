import unittest
from types import SimpleNamespace
from unittest.mock import patch

from primer_designer_app.utils.insilico_analysis import (
    insilico_reference_description,
    prepare_context_path,
    process_primer_pair,
)
from primer_designer_app.utils.primer_utils import (
    INSILICO_ERROR,
    INSILICO_OK,
    INSILICO_OK_EMPTY,
    PrimerPairResult,
)


def _pair():
    return PrimerPairResult(
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


class PrepareContextPathTests(unittest.TestCase):
    def test_grch37_uses_lift37_transcriptome(self):
        settings = SimpleNamespace(reference_genome="GRCh37")
        paths = prepare_context_path(settings)
        self.assertIn("gencode.v37lift37.transcripts.fa.gz", paths["cdna"])
        self.assertIn("Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz", paths["dna"])

    def test_grch38_uses_v49_transcriptome(self):
        settings = SimpleNamespace(reference_genome="GRCh38")
        paths = prepare_context_path(settings)
        self.assertIn("gencode.v49.transcripts.fa.gz", paths["cdna"])
        self.assertIn("Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", paths["dna"])

    def test_unsupported_genome_raises(self):
        settings = SimpleNamespace(reference_genome="hg18")
        with self.assertRaisesRegex(ValueError, "Unsupported reference genome"):
            prepare_context_path(settings)


class InsilicoReferenceDescriptionTests(unittest.TestCase):
    def test_disabled_returns_empty_string(self):
        settings = SimpleNamespace(
            do_insilico_pcr=False, reference_genome="GRCh37", context="genomic"
        )
        self.assertEqual(insilico_reference_description(settings), "")

    def test_genomic_context_mentions_dna_reference(self):
        settings = SimpleNamespace(
            do_insilico_pcr=True, reference_genome="GRCh37", context="genomic"
        )
        description = insilico_reference_description(settings)
        self.assertIn("Genomic amplicons", description)
        self.assertIn("Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz", description)

    def test_transcriptomic_context_mentions_cdna_reference(self):
        settings = SimpleNamespace(
            do_insilico_pcr=True, reference_genome="GRCh38", context="transcriptomic"
        )
        description = insilico_reference_description(settings)
        self.assertIn("Transcriptomic amplicons", description)
        self.assertIn("gencode.v49.transcripts.fa.gz", description)

    def test_unknown_context_returns_empty_string(self):
        settings = SimpleNamespace(
            do_insilico_pcr=True, reference_genome="GRCh37", context="bogus"
        )
        self.assertEqual(insilico_reference_description(settings), "")


class ProcessPrimerPairTests(unittest.TestCase):
    @patch("primer_designer_app.utils.insilico_analysis.run_dicey")
    def test_dicey_failure_returns_error_status(self, mock_run_dicey, *_):
        mock_run_dicey.return_value = None
        with self._tmp_dir() as temp_dir:
            result = process_primer_pair(_pair(), 0, "ref.fa.gz", temp_dir)
        self.assertEqual(result["insilico_status"], INSILICO_ERROR)
        self.assertEqual(result["amplicons"], [])
        self.assertIsNotNone(result["insilico_error_detail"])

    @patch("primer_designer_app.utils.insilico_analysis.run_dicey")
    def test_no_amplicons_found_returns_ok_empty(self, mock_run_dicey):
        mock_run_dicey.return_value = {"data": {"amplicons": []}}
        with self._tmp_dir() as temp_dir:
            result = process_primer_pair(_pair(), 0, "ref.fa.gz", temp_dir)
        self.assertEqual(result["insilico_status"], INSILICO_OK_EMPTY)
        self.assertEqual(result["amplicons"], [])
        self.assertEqual(result["insilico_seq"], "")

    @patch("primer_designer_app.utils.insilico_analysis.run_dicey")
    def test_single_matching_amplicon_sets_insilico_seq(self, mock_run_dicey):
        mock_run_dicey.return_value = {
            "data": {
                "amplicons": [
                    {"ForPos": 10, "RevEnd": 99, "Seq": "ACGTACGT"},
                ]
            }
        }
        pair = _pair()
        with self._tmp_dir() as temp_dir:
            result = process_primer_pair(pair, 0, "ref.fa.gz", temp_dir)
        self.assertEqual(result["insilico_status"], INSILICO_OK)
        self.assertEqual(len(result["amplicons"]), 1)
        self.assertEqual(result["insilico_seq"], "ACGTACGT")

    @patch("primer_designer_app.utils.insilico_analysis.run_dicey")
    def test_multiple_amplicons_returns_ok_without_insilico_seq(self, mock_run_dicey):
        mock_run_dicey.return_value = {
            "data": {
                "amplicons": [
                    {"ForPos": 10, "RevEnd": 99, "Seq": "ACGT"},
                    {"ForPos": 500, "RevEnd": 600, "Seq": "TTTT"},
                ]
            }
        }
        with self._tmp_dir() as temp_dir:
            result = process_primer_pair(_pair(), 0, "ref.fa.gz", temp_dir)
        self.assertEqual(result["insilico_status"], INSILICO_OK)
        self.assertEqual(len(result["amplicons"]), 2)
        self.assertEqual(result["insilico_seq"], "")

    def _tmp_dir(self):
        import tempfile

        return tempfile.TemporaryDirectory()


if __name__ == "__main__":
    unittest.main()
