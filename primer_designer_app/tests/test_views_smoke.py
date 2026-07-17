from django.test import Client, TestCase
from django.urls import reverse

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
