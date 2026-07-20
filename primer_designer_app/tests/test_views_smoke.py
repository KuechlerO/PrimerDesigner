from django.test import Client, TestCase
from django.urls import reverse
from unittest.mock import patch

from primer_designer_app.utils.primer_utils import INSILICO_OK, PrimerPairResult
from primer_designer_app.views.view_utils import (
    DOCX_CONTENT_TYPE,
    download_docx_report,
)


class StructuralVariantIndexViewTests(TestCase):
    def setUp(self):
        self.client = Client()

    def test_get_renders_index_page(self):
        response = self.client.get(
            reverse("primer_designer_app:structural_variants_index")
        )
        self.assertEqual(response.status_code, 200)

    def test_post_with_invalid_start_position_returns_400_error_page(self):
        response = self.client.post(
            reverse("primer_designer_app:structural_variants_index"),
            data={
                "sv_chromosome": "1",
                "sv_start_position": "not-a-number",
                "sv_end_position": "5000",
                "sv_type": "deletion",
            },
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Invalid input", response.content)
        self.assertIn(b"valid integer", response.content)


class AlleleSpecificPrimersOverviewViewTests(TestCase):
    def setUp(self):
        self.client = Client()

    def test_get_returns_method_not_allowed(self):
        response = self.client.get(
            reverse("primer_designer_app:allele_specific_primers_overview")
        )
        self.assertEqual(response.status_code, 405)
        self.assertIn(b"Method not allowed", response.content)


class SnvIndelPrimersOverviewViewTests(TestCase):
    def setUp(self):
        self.client = Client()

    def test_post_with_no_recognizable_fields_returns_400(self):
        response = self.client.post(
            reverse("primer_designer_app:snv_indel_primers_overview"), data={}
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Invalid input", response.content)
        self.assertIn(b"No recognizable input field found", response.content)


class DownloadDocxReportHelperTests(TestCase):
    def test_sets_content_disposition_and_docx_content_type(self):
        buffer = b"fake docx bytes"
        response = download_docx_report(buffer, "my_report.docx")

        self.assertEqual(
            response["Content-Disposition"], "attachment; filename=my_report.docx"
        )
        self.assertEqual(response["Content-Type"], DOCX_CONTENT_TYPE)
        self.assertEqual(response.content, buffer)


class SilicoPcrIndexViewTests(TestCase):
    def setUp(self):
        self.client = Client()
        self.url = reverse("primer_designer_app:silico_pcr_index")

    def test_get_renders_index_page(self):
        response = self.client.get(self.url)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "SilicoPCR")
        self.assertContains(response, "forward_primer")

    def test_post_missing_primer_returns_400(self):
        response = self.client.post(
            self.url,
            data={
                "forward_primer": "",
                "reverse_primer": "ATGCGATCGATCGATCGATC",
                "reference-genome": "GRCh37",
                "amplicon-check": "genome",
            },
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Invalid input", response.content)
        self.assertIn(b"Forward primer", response.content)

    def test_post_with_mocked_dicey_shows_results(self):
        pair = PrimerPairResult(
            index=0,
            left_seq="ATGCGATCGATCGATCGATC",
            right_seq="TACGTAGCTAGCTAGCTAGC",
            penalty=0.0,
            product_size=0,
            amplicons=[{"Length": 200, "Chrom": "1", "ForPos": 100, "RevEnd": 300}],
            insilico_status=INSILICO_OK,
        )

        def fake_run(**kwargs):
            from types import SimpleNamespace

            settings = SimpleNamespace(
                reference_genome="GRCh37",
                context="genomic",
                do_insilico_pcr=True,
            )
            return pair, settings, "Genomic amplicons were identified."

        with patch(
            "primer_designer_app.views.silico_pcr.run_silico_pcr",
            side_effect=fake_run,
        ):
            response = self.client.post(
                self.url,
                data={
                    "forward_primer": "ATGCGATCGATCGATCGATC",
                    "reverse_primer": "TACGTAGCTAGCTAGCTAGC",
                    "reference-genome": "GRCh37",
                    "amplicon-check": "genome",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "SilicoPCR Results")
        self.assertContains(response, "Show details")
        self.assertContains(response, "ATGCGATCGATCGATCGATC")
