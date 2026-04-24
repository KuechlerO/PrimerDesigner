from django.shortcuts import render
from primer_designer_app.views.view_utils import build_primer_settings
from primer_designer_app.utils.sv_utils import (
    build_structural_variant_info_from_request,
    design_structural_variant_primers,
)


def _build_template_form_data(request) -> dict:
    return {
        "reference_genome": request.POST.get("reference-genome", "GRCh37"),
        "usecase": request.POST.get("usecase", "qPCR"),
        "sv_chromosome": request.POST.get("sv_chromosome", ""),
        "sv_start_position": request.POST.get("sv_start_position", ""),
        "sv_end_position": request.POST.get("sv_end_position", ""),
        "sv_type": request.POST.get("sv_type", ""),
        "tm": request.POST.get("tm", "60"),
        "gc_content": request.POST.get("gc_content", "50"),
        "product_size_min": request.POST.get("product_size_min", "100"),
        "product_size_max": request.POST.get("product_size_max", "300"),
        "max_poly_X": request.POST.get("max_poly_X", "4"),
    }


def index(request):
    context = {
        "is_structural_variant": True,
    }

    if request.method == "POST":
        context["form_data"] = _build_template_form_data(request)

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
