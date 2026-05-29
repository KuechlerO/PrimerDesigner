# primer_designer_app/middleware.py
from django.shortcuts import render
from requests import HTTPError, RequestException

from .exceptions import (
    PrimerDesignerError,
    InvalidTranscriptIdError,
    InvalidTranscriptVersionError,
    InvalidTranscriptInputError,
    ExonExonJunctionError,
    InvalidReferenceSequenceError,
    NoPrimerPairsFoundError,
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
                "primer_designer_app/error_handling.html",
                {
                    "title": "Invalid Transcript ID",
                    "message": str(exception),
                    "hint": "Please enter a valid Ensembl transcript ID starting with ENST.",
                    "status_code": 400,
                },
                status=400,
            )

        if isinstance(exception, InvalidTranscriptVersionError):
            return render(
                request,
                "primer_designer_app/error_handling.html",
                {
                    "title": "Transcript version mismatch",
                    "message": str(exception),
                    "hint": "Please check whether the transcript version is correct.",
                    "status_code": 400,
                },
                status=400,
            )

        if isinstance(exception, InvalidTranscriptInputError):
            return render(
                request,
                "primer_designer_app/error_handling.html",
                {
                    "title": "Invalid transcript input",
                    "message": str(exception),
                    "hint": "Please provide a complete transcript-based input.",
                    "status_code": 400,
                },
                status=400,
            )

        if isinstance(exception, ExonExonJunctionError):
            return render(
                request,
                "primer_designer_app/error_handling.html",
                {
                    "title": "Unsupported variant location",
                    "message": str(exception),
                    "hint": "Variants affecting exon-exon junctions are currently not supported.",
                    "status_code": 400,
                },
                status=400,
            )

        if isinstance(exception, InvalidReferenceSequenceError):
            return render(
                request,
                "primer_designer_app/error_handling.html",
                {
                    "title": "Reference sequence not usable",
                    "message": str(exception),
                    "hint": (
                        "Check chromosome, position, and genome build (GRCh37 vs GRCh38). "
                        "For genomic SNV input, also provide the alternate base in the "
                        '"New base" field (e.g. Chr13:2655000 with new base A).'
                    ),
                    "status_code": 400,
                },
                status=400,
            )

        if isinstance(exception, NoPrimerPairsFoundError):
            return render(
                request,
                "primer_designer_app/error_handling.html",
                {
                    "title": "No primer pairs found",
                    "message": str(exception),
                    "hint": (
                        "Relax primer constraints (length, Tm, GC), widen the product size "
                        "range, or increase target padding in the primer settings."
                    ),
                    "status_code": 400,
                },
                status=400,
            )

        if isinstance(exception, PrimerDesignerError):
            return render(
                request,
                "primer_designer_app/error_handling.html",
                {
                    "title": "Primer design error",
                    "message": str(exception),
                    "hint": "Please review your input and primer settings, then try again.",
                    "status_code": 400,
                },
                status=400,
            )

        if isinstance(exception, HTTPError):
            return render(
                request,
                "primer_designer_app/error_handling.html",
                {
                    "title": "External service error",
                    "message": str(exception),
                    "hint": "The Ensembl service returned an HTTP error. Please try again later.",
                    "status_code": 503,
                },
                status=503,
            )

        if isinstance(exception, RequestException):
            return render(
                request,
                "primer_designer_app/error_handling.html",
                {
                    "title": "Connection error",
                    "message": str(exception),
                    "hint": "The external service could not be reached. Please try again later.",
                    "status_code": 503,
                },
                status=503,
            )

        return None
