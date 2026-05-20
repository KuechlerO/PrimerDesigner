from django.shortcuts import render
from primer_designer_app.views.view_utils import (
    build_form_data_from_request,
    build_primer_settings,
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

            context["structural_variant_info"] = structural_variant_info
            context["sv_results"] = sv_results

        except Exception as exc:
            context["error_message"] = str(exc)

    return render(
        request,
        "primer_designer_app/structural_variant.html",
        context,
    )
