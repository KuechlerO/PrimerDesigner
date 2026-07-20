from django.shortcuts import render

from primer_designer_app.exceptions import InvalidInputError
from primer_designer_app.models import DesignResultsSummary
from primer_designer_app.utils.doc_utils import create_structural_variant_primer_report
from primer_designer_app.views.view_utils import (
    build_form_data_from_request,
    build_primer_settings,
    download_docx_report,
)
from primer_designer_app.utils.sv_utils import (
    build_structural_variant_info_from_request,
    design_structural_variant_primers,
)


def index(request):
    sv_fields = {}
    if request.method == "POST":
        sv_fields = {
            "sv_chromosome": request.POST.get("sv_chromosome", ""),
            "sv_start_position": request.POST.get("sv_start_position", ""),
            "sv_end_position": request.POST.get("sv_end_position", ""),
            "sv_type": request.POST.get("sv_type", ""),
        }

    context = {
        "is_structural_variant": True,
        "form_data": build_form_data_from_request(request, **sv_fields),
    }

    if request.method == "POST":
        try:
            primer_settings = build_primer_settings(request)
            primer_settings.do_insilico_pcr = False

            structural_variant_info = build_structural_variant_info_from_request(
                request
            )

            sv_results = design_structural_variant_primers(
                structural_variant_info=structural_variant_info,
                primer_settings=primer_settings,
            )

            result_summary = DesignResultsSummary()
            result_summary.save_structural_variant_results(
                primer_settings, structural_variant_info, sv_results
            )

            context["structural_variant_info"] = structural_variant_info
            context["sv_results"] = sv_results
            context["report_uuid"] = result_summary.id
        except ValueError as exc:
            raise InvalidInputError(str(exc)) from exc

    return render(
        request,
        "primer_designer_app/structural_variant.html",
        context,
    )


def generate_report(request, uuid):
    """Download a Word report with all SV primer pairs (no sequence figure)."""
    design_summary = DesignResultsSummary.objects.get(id=uuid)
    if not design_summary.is_structural_variant_design():
        raise InvalidInputError("Not a structural variant design result.")

    doc_buffer = create_structural_variant_primer_report(design_summary)
    return download_docx_report(doc_buffer, f"sv_primer_report_{uuid}.docx")
