from django.test import RequestFactory, TestCase

from primer_designer_app.exceptions import (
    InvalidInputError,
    InvalidTranscriptIdError,
    NoPrimerPairsFoundError,
)
from primer_designer_app.middleware import PrimerDesignerErrorMiddleware


class PrimerDesignerErrorMiddlewareTests(TestCase):
    def setUp(self):
        self.factory = RequestFactory()
        self.middleware = PrimerDesignerErrorMiddleware(
            get_response=lambda request: None
        )

    def _request(self):
        return self.factory.get("/some-url/")

    def test_invalid_input_error_returns_400_with_expected_title(self):
        response = self.middleware.process_exception(
            self._request(), InvalidInputError("bad field")
        )
        self.assertIsNotNone(response)
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Invalid input", response.content)
        self.assertIn(b"bad field", response.content)

    def test_no_primer_pairs_found_error_returns_400(self):
        response = self.middleware.process_exception(
            self._request(),
            NoPrimerPairsFoundError("Primer3 did not return any primer pairs."),
        )
        self.assertIsNotNone(response)
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"No primer pairs found", response.content)

    def test_invalid_transcript_id_error_returns_400(self):
        response = self.middleware.process_exception(
            self._request(), InvalidTranscriptIdError("ENST is not valid")
        )
        self.assertIsNotNone(response)
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"Invalid Transcript ID", response.content)

    def test_unhandled_exception_type_returns_none(self):
        response = self.middleware.process_exception(self._request(), KeyError("oops"))
        self.assertIsNone(response)

    def test_call_delegates_to_get_response(self):
        sentinel = object()
        middleware = PrimerDesignerErrorMiddleware(
            get_response=lambda request: sentinel
        )
        self.assertIs(middleware(self._request()), sentinel)
