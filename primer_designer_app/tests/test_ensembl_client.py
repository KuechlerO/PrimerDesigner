from unittest.mock import MagicMock, patch

import requests
from django.test import SimpleTestCase

from primer_designer_app.utils.ensembl_client import (
    TIMEOUT_SEQUENCE,
    EnsemblClient,
    server37,
    server38,
)


class EnsemblClientConfigTests(SimpleTestCase):
    def test_grch38_uses_https_rest_server(self):
        client = EnsemblClient(ref_genome="GRCh38")
        self.assertEqual(client.server, "https://rest.ensembl.org")
        self.assertEqual(server38, "https://rest.ensembl.org")

    def test_grch37_uses_archive_https_server(self):
        client = EnsemblClient(ref_genome="GRCh37")
        self.assertEqual(client.server, "https://grch37.rest.ensembl.org")
        self.assertEqual(server37, "https://grch37.rest.ensembl.org")

    def test_retry_adapter_mounted_for_http_and_https(self):
        client = EnsemblClient(ref_genome="GRCh38")
        self.assertIn("https://", client.session.adapters)
        self.assertIn("http://", client.session.adapters)
        https_adapter = client.session.adapters["https://"]
        http_adapter = client.session.adapters["http://"]
        self.assertIsNotNone(https_adapter.max_retries)
        self.assertEqual(https_adapter.max_retries.total, 2)
        self.assertEqual(http_adapter.max_retries.total, 2)
        self.assertNotIn(500, https_adapter.max_retries.status_forcelist)

    def test_split_transcript_id_strips_version(self):
        client = EnsemblClient(ref_genome="GRCh38")
        base_id, version = client.split_transcript_id("ENST00000424479.7")
        self.assertEqual(base_id, "ENST00000424479")
        self.assertEqual(version, 7)

    def test_split_transcript_id_without_version(self):
        client = EnsemblClient(ref_genome="GRCh38")
        base_id, version = client.split_transcript_id("ENST00000424479")
        self.assertEqual(base_id, "ENST00000424479")
        self.assertIsNone(version)

    @patch.object(EnsemblClient, "__init__", lambda self, ref_genome="GRCh38": None)
    def test_grch37_prefers_current_ensembl_gene_symbol(self):
        client = EnsemblClient()
        client.server = "https://grch37.rest.ensembl.org"
        client.session = MagicMock()

        lookup_transcript = MagicMock()
        lookup_transcript.raise_for_status = MagicMock()
        lookup_transcript.json.return_value = {
            "version": 2,
            "Parent": "ENSG00000147535",
            "display_name": "PPAPDC1B-001",
        }

        lookup_gene_37 = MagicMock()
        lookup_gene_37.raise_for_status = MagicMock()
        lookup_gene_37.json.return_value = {"display_name": "PPAPDC1B"}

        lookup_gene_38 = MagicMock()
        lookup_gene_38.raise_for_status = MagicMock()
        lookup_gene_38.json.return_value = {"display_name": "PLPP5"}

        def get_side_effect(url, **kwargs):
            if "grch37" in url and "ENSG" in url:
                return lookup_gene_37
            if "rest.ensembl.org/lookup/id/ENSG" in url and "grch37" not in url:
                return lookup_gene_38
            return lookup_transcript

        client.session.get.side_effect = get_side_effect

        symbol, gene_id, full_id = client.get_gene_symbol_for_transcriptID(
            "ENST00000424479.2"
        )
        self.assertEqual(symbol, "PLPP5")
        self.assertEqual(gene_id, "ENSG00000147535")
        self.assertEqual(full_id, "ENST00000424479.2")

    @patch.object(EnsemblClient, "__init__", lambda self, ref_genome="GRCh38": None)
    def test_gene_lookup_500_falls_back_to_transcript_display_name(self):
        client = EnsemblClient()
        client.server = "https://rest.ensembl.org"
        client.session = MagicMock()

        lookup_transcript = MagicMock()
        lookup_transcript.raise_for_status = MagicMock()
        lookup_transcript.json.return_value = {
            "version": 7,
            "Parent": "ENSG00000147535",
            "display_name": "PLPP5-203",
        }

        def get_side_effect(url, **kwargs):
            if "ENSG" in url:
                raise requests.HTTPError(
                    "500 Server Error: Internal Server Error for url: "
                    "https://rest.ensembl.org/lookup/id/ENSG00000147535"
                )
            return lookup_transcript

        client.session.get.side_effect = get_side_effect
        symbol, gene_id, full_id = client.get_gene_symbol_for_transcriptID(
            "ENST00000424479.7"
        )
        self.assertEqual(symbol, "PLPP5")
        self.assertEqual(gene_id, "ENSG00000147535")
        self.assertEqual(full_id, "ENST00000424479.7")

    @patch.object(EnsemblClient, "__init__", lambda self, ref_genome="GRCh38": None)
    def test_sequence_request_uses_base_id_and_sequence_timeout(self):
        client = EnsemblClient()
        client.server = "https://rest.ensembl.org"
        client.session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "ACGT"
        mock_response.raise_for_status = MagicMock()
        client.session.get.return_value = mock_response

        seq = client.get_transcript_sequence("ENST00000424479.7", "cdna")

        self.assertEqual(seq, "ACGT")
        args, kwargs = client.session.get.call_args
        self.assertIn("/sequence/id/ENST00000424479?", args[0])
        self.assertNotIn("mask_feature", args[0])
        self.assertNotIn(".7", args[0].split("/sequence/id/")[1].split("?")[0])
        self.assertEqual(kwargs["timeout"], TIMEOUT_SEQUENCE)

    @patch.object(EnsemblClient, "__init__", lambda self, ref_genome="GRCh38": None)
    def test_sequence_retries_after_transient_500(self):
        client = EnsemblClient()
        client.server = "https://rest.ensembl.org"
        client.session = MagicMock()

        fail = MagicMock()
        fail.status_code = 500
        fail.text = "<!doctype html>error"
        fail.raise_for_status.side_effect = requests.HTTPError("500")

        ok = MagicMock()
        ok.status_code = 200
        ok.text = "ACGTACGT"
        ok.raise_for_status = MagicMock()

        # plain fail, json fail, then after delay plain ok
        client.session.get.side_effect = [fail, fail, ok]

        with patch("primer_designer_app.utils.ensembl_client.time.sleep"):
            with patch(
                "primer_designer_app.utils.ensembl_client._SEQUENCE_RETRY_DELAYS_SEC",
                (0.0, 0.0),
            ):
                seq = client.get_transcript_sequence("ENST00000424479.7", "cdna")

        self.assertEqual(seq, "ACGTACGT")
        self.assertGreaterEqual(client.session.get.call_count, 3)
