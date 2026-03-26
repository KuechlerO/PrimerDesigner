# primer_designer_app/middleware.py
from django.shortcuts import render
from requests import HTTPError, RequestException

from .exceptions import (
    PrimerDesignerError,
    InvalidTranscriptIdError,
    InvalidTranscriptVersionError,
    InvalidTranscriptInputError,
    ExonExonJunctionError,
)

class PrimerDesignerErrorMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        return self.get_response(request)

    def process_exception(self, request, exception):
        if isinstance(exception, InvalidTranscriptIdError):
            return render(
                request,
                'primer_designer_app/error_handling.html',
                {
                    'title': 'Invalid Transcript ID',
                    'message': str(exception),
                    'hint': 'Please enter a valid Ensembl transcript ID starting with ENST.',
                    'status_code': 400,
                },
                status=400,
            )

        if isinstance(exception, InvalidTranscriptVersionError):
            return render(
                request,
                'primer_designer_app/error_handling.html',
                {
                    'title': 'Transcript version mismatch',
                    'message': str(exception),
                    'hint': 'Please check whether the transcript version is correct.',
                    'status_code': 400,
                },
                status=400,
            )

        if isinstance(exception, InvalidTranscriptInputError):
            return render(
                request,
                'primer_designer_app/error_handling.html',
                {
                    'title': 'Invalid transcript input',
                    'message': str(exception),
                    'hint': 'Please provide a complete transcript-based input.',
                    'status_code': 400,
                },
                status=400,
            )

        if isinstance(exception, ExonExonJunctionError):
            return render(
                request,
                'primer_designer_app/error_handling.html',
                {
                    'title': 'Unsupported variant location',
                    'message': str(exception),
                    'hint': 'Variants affecting exon-exon junctions are currently not supported.',
                    'status_code': 400,
                },
                status=400,
            )

        if isinstance(exception, HTTPError):
            return render(
                request,
                'primer_designer_app/error_handling.html',
                {
                    'title': 'External service error',
                    'message': str(exception),
                    'hint': 'The Ensembl service returned an HTTP error. Please try again later.',
                    'status_code': 503,
                },
                status=503,
            )

        if isinstance(exception, RequestException):
            return render(
                request,
                'primer_designer_app/error_handling.html',
                {
                    'title': 'Connection error',
                    'message': str(exception),
                    'hint': 'The external service could not be reached. Please try again later.',
                    'status_code': 503,
                },
                status=503,
            )

        return None
