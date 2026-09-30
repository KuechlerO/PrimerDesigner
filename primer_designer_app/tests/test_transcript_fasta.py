import os
import subprocess
import tempfile
import unittest

from django.test import override_settings

from primer_designer_app.exceptions import InvalidReferenceSequenceError
from primer_designer_app.utils.transcript_fasta import (
    fetch_transcript_sequence,
    transcript_fasta_path,
)


class TranscriptFastaPathTests(unittest.TestCase):
    @override_settings(REFERENCE_DATA_DIR="/tmp/refs")
    def test_resolves_grch38_cdna_and_cds(self):
        self.assertEqual(
            transcript_fasta_path("GRCh38", "cdna"),
            "/tmp/refs/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        )
        self.assertEqual(
            transcript_fasta_path("GRCh38", "cds"),
            "/tmp/refs/Homo_sapiens.GRCh38.cds.all.fa.gz",
        )

    @override_settings(REFERENCE_DATA_DIR="/tmp/refs")
    def test_resolves_grch37(self):
        self.assertEqual(
            transcript_fasta_path("GRCh37", "cdna"),
            "/tmp/refs/Homo_sapiens.GRCh37.cdna.all.fa.gz",
        )

    def test_rejects_unsupported_seq_type(self):
        with self.assertRaises(ValueError):
            transcript_fasta_path("GRCh38", "protein")


class TranscriptFastaFetchTests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.ref_dir = self._tmpdir.name
        self.fasta_name = "Homo_sapiens.GRCh38.cdna.all.fa.gz"
        self.fasta_path = os.path.join(self.ref_dir, self.fasta_name)

        plain = os.path.join(self.ref_dir, "tmp.fa")
        with open(plain, "w", encoding="utf-8") as fh:
            fh.write(">ENST00000424479.7\n")
            fh.write("ATGCGATCGATCG\n")
            fh.write(">ENST00000000001.1\n")
            fh.write("AAAA\n")

        # Prefer bgzip if available; gzip is fine for samtools faidx on many builds.
        bgzip = subprocess.run(
            ["bgzip", "-c", plain],
            check=False,
            capture_output=True,
        )
        if bgzip.returncode == 0 and bgzip.stdout:
            with open(self.fasta_path, "wb") as out:
                out.write(bgzip.stdout)
        else:
            import gzip

            with open(plain, "rb") as src, gzip.open(self.fasta_path, "wb") as out:
                out.write(src.read())

        faidx = subprocess.run(
            ["samtools", "faidx", self.fasta_path],
            check=False,
            capture_output=True,
            text=True,
        )
        if faidx.returncode != 0:
            self.skipTest(f"samtools faidx unavailable: {faidx.stderr}")

    def tearDown(self):
        self._tmpdir.cleanup()

    def test_exact_id_hit(self):
        with override_settings(REFERENCE_DATA_DIR=self.ref_dir):
            seq = fetch_transcript_sequence(
                "ENST00000424479.7", "cdna", "GRCh38"
            )
        self.assertEqual(seq, "ATGCGATCGATCG")

    def test_missing_id(self):
        with override_settings(REFERENCE_DATA_DIR=self.ref_dir):
            with self.assertRaisesRegex(
                InvalidReferenceSequenceError, "not found in local CDNA FASTA"
            ):
                fetch_transcript_sequence(
                    "ENST00000424479.2", "cdna", "GRCh38"
                )

    def test_missing_file(self):
        empty = tempfile.mkdtemp()
        try:
            with override_settings(REFERENCE_DATA_DIR=empty):
                with self.assertRaisesRegex(
                    InvalidReferenceSequenceError, "reference file not found"
                ):
                    fetch_transcript_sequence(
                        "ENST00000424479.7", "cdna", "GRCh38"
                    )
        finally:
            os.rmdir(empty)

    def test_missing_fai(self):
        os.remove(self.fasta_path + ".fai")
        with override_settings(REFERENCE_DATA_DIR=self.ref_dir):
            with self.assertRaisesRegex(
                InvalidReferenceSequenceError, "Missing FASTA index"
            ):
                fetch_transcript_sequence(
                    "ENST00000424479.7", "cdna", "GRCh38"
                )


if __name__ == "__main__":
    unittest.main()
