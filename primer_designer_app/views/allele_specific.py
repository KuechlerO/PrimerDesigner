from django.shortcuts import render
from django.http import HttpResponse

import logging

from primer_designer_app.models import DesignResultsSummary
from primer_designer_app.utils.doc_utils import create_primer_report
from primer_designer_app.utils.primer_utils import PrimerSearchResults
from primer_designer_app.views.view_utils import (
    _get_post,
    build_form_data_from_request,
    build_primer_settings,
    handle_allele_specific_input,
)

logger = logging.getLogger(__name__)


def index(request):
    """Index page for allele-specific PCR mode (AS-PCR / ARMS-PCR)."""
    return render(
        request,
        "primer_designer_app/allele_specific_index.html",
        {
            "form_data": build_form_data_from_request(request),
            "is_allele_specific": True,
        },
    )


def primers_overview(request, uuid=None):
    """
    AS-PCR overview.

    Current implementation persists WT+MUT designs, but reuses the existing
    SNV/Indel results template as an interim UI (showing WT set only).
    """
    if request.method != "POST":
        return HttpResponse("Method not allowed", status=405)

    primer_settings_obj = build_primer_settings(request)
    primer_settings_obj.do_insilico_pcr = False
    primer_settings_obj.check_known_snps = False

    try:
        new_uuid = handle_allele_specific_input(request, primer_settings_obj)
    except Exception as exc:
        return HttpResponse(str(exc), status=400)

    designResults_obj = DesignResultsSummary.objects.get(id=new_uuid)
    var_info = designResults_obj.get_variant_info()
    as_res = designResults_obj.get_primer_search_results()
    if not isinstance(as_res, dict) or as_res.get("design_type") != "allele_specific":
        return HttpResponse("AS-PCR design data missing or invalid.", status=500)

    wt_results = as_res["wt"]
    mut_results = as_res["mut"]

    if not (wt_results.primer_pairs and mut_results.primer_pairs):
        return HttpResponse(
            "AS-PCR: Primer3 did not return any primer pairs for WT and/or MUT.",
            status=400,
        )

    from primer_designer_app.utils.helpers import (
        html_visualize_sequence_allele_specific,
        create_hgvs_notation,
        vcf_hits_json_for_display,
    )

    highlighted_seq_wt, display_offset_wt, display_length_wt, wt_display_chunks = (
        html_visualize_sequence_allele_specific(
            designResults_obj.primer_settings,
            var_info,
            wt_results.primer_pairs[0],
            allele="wt",
            all_primer_pairs=wt_results.primer_pairs,
        )
    )
    highlighted_seq_mut, display_offset_mut, display_length_mut, mut_display_chunks = (
        html_visualize_sequence_allele_specific(
            designResults_obj.primer_settings,
            var_info,
            mut_results.primer_pairs[0],
            allele="mut",
            all_primer_pairs=mut_results.primer_pairs,
        )
    )

    hgvs_info = create_hgvs_notation(var_info)

    vcf_applied = (designResults_obj.variant_info_data or {}).get(
        "vcf_applied_variants"
    ) or []
    vcf_hits_json_wt = vcf_hits_json_for_display(
        vcf_applied, display_offset_wt, display_length_wt
    )
    vcf_hits_json_mut = vcf_hits_json_for_display(
        vcf_applied, display_offset_mut, display_length_mut
    )

    return render(
        request,
        "primer_designer_app/allele_specific_results.html",
        {
            "result_sum_obj": designResults_obj,
            "hgvs_info": hgvs_info,
            "common_reverse_primer": as_res.get("common_reverse_primer") or "",
            "wt_results": wt_results,
            "mut_results": mut_results,
            "highlighted_sequence_wt": highlighted_seq_wt,
            "highlighted_sequence_mut": highlighted_seq_mut,
            "wt_display_chunks": wt_display_chunks,
            "mut_display_chunks": mut_display_chunks,
            "vcf_hits_json_wt": vcf_hits_json_wt,
            "vcf_hits_json_mut": vcf_hits_json_mut,
            "sequence_display_offset_wt": display_offset_wt,
            "sequence_display_offset_mut": display_offset_mut,
        },
    )


def primer_details(request, uuid):
    return HttpResponse(
        "Allele-specific PCR details view not implemented yet.", status=501
    )


def generate_report(request, uuid, selected_primer_index: int):
    """
    View to generate the final report after the user has decided on a primer pair.
    """
    designResults_obj = DesignResultsSummary.objects.get(id=uuid)
    doc_buffer = create_primer_report(designResults_obj, selected_primer_index)
    response = HttpResponse(
        doc_buffer,
        content_type="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
    )
    response["Content-Disposition"] = (
        f"attachment; filename=report_{uuid}_{selected_primer_index}.docx"
    )
    return response
