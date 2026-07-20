from django.shortcuts import render

from primer_designer_app.views.view_utils import build_form_data_from_request
from primer_designer_app.utils.silico_pcr_utils import run_silico_pcr


def index(request):
    """SilicoPCR: paste forward/reverse primers and run Dicey amplicon search."""
    primer_fields = {}
    if request.method == "POST":
        primer_fields = {
            "forward_primer": request.POST.get("forward_primer", ""),
            "reverse_primer": request.POST.get("reverse_primer", ""),
        }

    form_data = build_form_data_from_request(request, **primer_fields)
    # SilicoPCR always runs Dicey; default context is genome (not "none").
    if not request.POST.get("amplicon-check"):
        form_data["amplicon_check"] = "genome"
    elif form_data.get("amplicon_check") == "none":
        form_data["amplicon_check"] = "genome"

    context = {
        "is_silico_pcr": True,
        "form_data": form_data,
    }

    if request.method == "POST":
        pair, settings, reference_note = run_silico_pcr(
            forward_raw=request.POST.get("forward_primer", ""),
            reverse_raw=request.POST.get("reverse_primer", ""),
            reference_genome=request.POST.get("reference-genome", "GRCh37"),
            amplicon_check=request.POST.get("amplicon-check", "genome"),
        )
        context["primer_pair"] = pair
        context["silico_settings"] = settings
        context["insilico_reference_note"] = reference_note
        context["show_results"] = True

    return render(request, "primer_designer_app/silico_pcr.html", context)
