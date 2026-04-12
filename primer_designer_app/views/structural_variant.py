from django.http import HttpResponse
from django.shortcuts import render, redirect



# def index(request):
#     # Return placeholder: Not yet implemented
#     return render(request, "primer_designer_app/structural_variant.html")



from primer_designer_app.views.view_utils import build_primer_settings
from primer_designer_app.utils.primer3_post import parse_primer3_overrides_from_post
from primer_designer_app.utils.sv_utils import (
    build_structural_variant_info_from_request,
    design_structural_variant_primers,
)


def _build_template_form_data(request) -> dict:
    return {
        'reference_genome': request.POST.get('reference-genome', 'GRCh37'),
        'sv_chromosome': request.POST.get('sv_chromosome', ''),
        'sv_start_position': request.POST.get('sv_start_position', ''),
        'sv_end_position': request.POST.get('sv_end_position', ''),
        'sv_type': request.POST.get('sv_type', ''),
        'tm': request.POST.get('tm', '60'),
        'gc_content': request.POST.get('gc_content', '50'),
        'product_size_min': request.POST.get('product_size_min', '100'),
        'product_size_max': request.POST.get('product_size_max', '300'),
        'max_poly_X': request.POST.get('max_poly_X', '4'),
    }


def index(request):
    context = {
        'form_data': {
            'reference_genome': 'GRCh37',
            'sv_chromosome': '',
            'sv_start_position': '',
            'sv_end_position': '',
            'sv_type': 'deletion',
            'tm': '60',
            'gc_content': '50',
            'product_size_min': '100',
            'product_size_max': '300',
            'max_poly_X': '4',
        }
    }

    if request.method == 'POST':
        context['form_data'] = _build_template_form_data(request)

        try:
            primer_settings = build_primer_settings(request)
            primer_settings.primer3_overrides = parse_primer3_overrides_from_post(request)
            primer_settings.do_insilico_pcr = False

            structural_variant_info = build_structural_variant_info_from_request(request)

            sv_results = design_structural_variant_primers(
                structural_variant_info=structural_variant_info,
                primer_settings=primer_settings,
            )

            context['structural_variant_info'] = structural_variant_info
            context['sv_results'] = sv_results

        except Exception as exc:
            context['error_message'] = str(exc)

    return render(
        request,
        'primer_designer_app/structural_variant.html',
        context,
    )
